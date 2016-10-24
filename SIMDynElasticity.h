// $Id$
//==============================================================================
//!
//! \file SIMDynElasticity.h
//!
//! \date Dec 04 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Dynamic simulation driver for elasticity problems with fracture.
//!
//==============================================================================

#ifndef _SIM_DYN_ELASTICITY_H_
#define _SIM_DYN_ELASTICITY_H_

#include "NewmarkSIM.h"
#include "SIMElasticity.h"
#include "InitialConditionHandler.h"
#include "FractureElasticityMonol.h"
#include "DataExporter.h"
#include "MultiPatchModelGenerator.h"
#include <fstream>


/*!
  \brief Driver class for dynamic elasticity problems with fracture.
*/

template<class Dim, class DynSIM=NewmarkSIM>
class SIMDynElasticity : public SIMElasticity<Dim>
{
public:
  //! \brief Default constructor.
  SIMDynElasticity(bool monolithic = false) : dSim(*this)
  {
    Dim::myHeading = "Elasticity solver";
    vtfStep = subIter = 0;
    if (monolithic)
    {
      ++Dim::nf[0]; // Account for the phase field as an unknown field variable
      phOrder = 2;  // Default 2nd order phase field when monolithic coupling
    }
    else
      phOrder = 0;
  }

  //! \brief Empty destructor.
  virtual ~SIMDynElasticity() {}

  //! \brief Prints out problem-specific data to the log stream.
  virtual void printProblem() const
  {
    static short int ncall = 0;
    if (++ncall == 1) // Avoiding infinite recursive calls
      dSim.printProblem();
    else
      this->SIMElasticity<Dim>::printProblem();
    --ncall;
  }

  //! \brief Registers fields for data output.
  void registerFields(DataExporter& exporter)
  {
    int flag = DataExporter::PRIMARY;
    if (!Dim::opt.pSolOnly)
      flag |= DataExporter::SECONDARY;
    exporter.registerField("u","solid displacement",DataExporter::SIM,flag);
    exporter.setFieldValue("u",this,&dSim.getSolution());
  }

  //! \brief Initializes the problem.
  bool init(const TimeStep&)
  {
    dSim.initPrm();
    dSim.initSol(3);
    if (phOrder > 1)
    {
      if (strcmp(DynSIM::inputContext,"nonlinearsolver"))
      {
        std::cerr <<" *** SIMDynElasticity::init: Monolithic dynamic simulation"
                  <<" is not available."<< std::endl;
        return false;
      }

      // Insert initial phase field solution 1.0 (undamaged material)
      size_t nndof = Dim::nf[0];
      Vector sol(this->getNoDOFs());
      for (size_t d = nndof-1; d < sol.size(); d += nndof)
        sol[d] = 1.0;
      dSim.setSolution(sol,0);
    }

    bool ok = this->setMode(SIM::INIT);
    this->setQuadratureRule(Dim::opt.nGauss[0],true);
    this->registerField("solution",dSim.getSolution());
    return SIM::setInitialConditions(*this) && ok;
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0) return true;

    return dSim.saveModel(geoBlk,nBlock,fileName);
  }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param nBlock Running result block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    double old = utl::zero_print_tol;
    utl::zero_print_tol = 1e-16;
    bool ok = this->savePoints(dSim.getSolution(),tp.time.t,tp.step);
    utl::zero_print_tol = old;

    if (!energFile.empty() && Dim::adm.getProcId() == 0 && tp.step > 0)
    {
      size_t i;
      Vector aLoad, react; // Integrate forces along specified boundary, if any
      this->getBoundaryForce(aLoad,dSim.getSolutions(),tp);
      this->getBoundaryReactions(react);

      // Write global energies and load resultants to the energyFile
      std::ofstream os(energFile, tp.step == 1 ? std::ios::out : std::ios::app);
      if (tp.step == 1)
      {
        os <<"#t eps_e external_energy eps+ eps- eps_b ";
        for (i = 0; i < aLoad.size(); i++) os <<" load_"<< char('X'+i);
        for (i = 0; i < react.size(); i++) os <<" react_"<< char('X'+i);
        os <<"\n";
      }
      os << std::setprecision(11) << std::setw(6) << std::scientific;
      os << tp.time.t;
      for (i = 0; i < gNorm.size(); i++) os <<" "<< gNorm[i];
      for (i = 0; i < aLoad.size(); i++) os <<" "<< utl::trunc(aLoad[i]);
      for (i = 0; i < react.size(); i++) os <<" "<< utl::trunc(react[i]);
      os <<"\n";
    }

    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0 || !ok)
      return ok;

    // Write primary and secondary (of requested) solution fields to VTF-file
    if (!dSim.saveStep(++vtfStep,nBlock,tp.time.t))
      return false;
    else if (tp.step < 1)
      return true;

    // Write projected solution fields to VTF-file
    if (!Dim::opt.project.empty())
      if (!this->writeGlvP(projSol,vtfStep,nBlock,110,
                           Dim::opt.project.begin()->second.c_str()))
        return false;

    // Write element norms to VTF-file
    return this->writeGlvN(eNorm,vtfStep,nBlock);
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp) { return dSim.advanceStep(tp,false); }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    if (Dim::msgLevel >= 1 && phOrder < 2)
      IFEM::cout <<"\n  Solving the elasto-dynamics problem...";

    if (dSim.solveStep(tp) != SIM::CONVERGED)
      return false;

    return this->postSolve(tp);
  }

  //! \brief Computes solution norms, etc. on the converged solution.
  bool postSolve(TimeStep& tp)
  {
    RealArray RF;
    if (this->getCurrentReactions(RF,dSim.getSolution()))
    {
      IFEM::cout <<"  Total reaction forces:          Sum(R) :";
      for (size_t i = 1; i < RF.size(); i++)
        IFEM::cout <<" "<< utl::trunc(RF[i]);
      double Ru = RF.front();
      if (utl::trunc(Ru) != 0.0)
        IFEM::cout <<"\n  displacement*reactions:          (R,u) : "<< Ru;
      IFEM::cout << std::endl;
    }

    // Update strain energy density for the converged solution
    this->setMode(SIM::RECOVERY);
    if (!this->assembleSystem(tp.time,dSim.getSolutions()))
      return false;

    // Project the secondary solution field onto the geometry basis
    if (!Dim::opt.project.empty())
      if (!this->project(projSol,dSim.getSolution(),
                         Dim::opt.project.begin()->first))
        return false;

    Vectors gNorms;
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    if (!this->solutionNorms(tp.time,dSim.getSolutions(),gNorms,&eNorm))
      return false;
    else if (!gNorms.empty())
    {
      gNorm = gNorms.front();
      if (gNorm.size() > 0 && utl::trunc(gNorm(1)) != 0.0)
        IFEM::cout <<"  Elastic strain energy:           eps_e : "
                   << gNorm(1) << std::endl;
      if (gNorm.size() > 4 && utl::trunc(gNorm(5)) != 0.0)
        IFEM::cout <<"  Bulk energy:                     eps_b : "
                   << gNorm(5)
                   <<"\n  Tensile & compressive energies         : "
                   << gNorm(3) <<" "<< gNorm(4) << std::endl;
      if (gNorm.size() > 5 && utl::trunc(gNorm(6)) != 0.0)
        IFEM::cout <<"  Dissipated energy:               eps_d : "
                   << gNorm(6) << std::endl;
      if (gNorm.size() > 1 && utl::trunc(gNorm(2)) != 0.0)
        IFEM::cout <<"  External energy: ((f,u^h)+(t,u^h))^0.5 : "
                   << (gNorm(2) < 0.0 ? -sqrt(-gNorm(2)) : sqrt(gNorm(2)))
                   << std::endl;
    }

    return true;
  }

  //! \brief Returns the tensile energy in gauss points.
  virtual const RealArray* getTensileEnergy() const
  {
    return static_cast<Elasticity*>(Dim::myProblem)->getTensileEnergy();
  }

  //! \brief Returns a const reference to the global norms.
  const Vector& getGlobalNorms() const { return gNorm; }

  //! \brief Parses sub-iteration parameters from an XML element.
  void parseSubiteration(const TiXmlElement* elem)
  {
    utl::getAttribute(elem,"type",subIter);
  }

  //! \brief Assigns the file name for global energy output (monolithic only).
  void setEnergyFile(const char* fName)
  {
    if (phOrder > 1 && fName)
    {
      energFile = fName;
      IFEM::cout <<"\tFile for global energy output: "<< energFile << std::endl;
    }
  }

  //! \brief Returns a const reference to current solution vector.
  const Vector& getSolution(int idx = 0) const { return dSim.getSolution(idx); }
  //! \brief Returns a const reference to the solution vectors.
  const Vectors& getSolutions() const { return dSim.getSolutions(); }

  //! \brief Updates the solution vectors.
  void setSolutions(const Vectors& dvec)
  {
    size_t nSol = dSim.getSolutions().size();
    for (size_t i = 0; i < nSol && i < dvec.size(); i++)
      dSim.setSolution(dvec[i],i);
  }

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] tp Time stepping parameters
  SIM::ConvStatus solveIteration(TimeStep& tp)
  {
    return subIter == 1 ? dSim.solveStep(tp) : dSim.solveIteration(tp);
  }

  //! \brief Returns the maximum number of iterations.
  int getMaxit() const { return dSim.getMaxit(); }

  //! \brief Checks whether an internal crack pressure has been specified.
  bool haveCrackPressure() const
  {
    FractureElasticity* fel = dynamic_cast<FractureElasticity*>(Dim::myProblem);
    return fel ? fel->getCrackPressure() != 0.0 : false;
  }

protected:
  ModelGenerator* createModelGenerator(const TiXmlElement* geo) const
  {
  IFEM::cout <<"  Using multi-patch model generator" << std::endl;
  if (Dim::dimension == 2)
    return new MultiPatchModelGenerator2D(geo);
  else
    return new MultiPatchModelGenerator3D(geo);
  }
  
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
    {
      if (phOrder > 1)
        Dim::myProblem = new FractureElasticityMonol(Dim::dimension,phOrder);
      else // Using the Voigt elasticity formulation by default
        Dim::myProblem = new FractureElasticityVoigt(Dim::dimension);
    }
    return static_cast<Elasticity*>(Dim::myProblem);
  }

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    bool result = true;
    static short int ncall = 0;
    if (++ncall == 1) // Avoiding infinite recursive calls
      result = dSim.parse(elem);
    else if (!strcasecmp(elem->Value(),"cahnhilliard") && phOrder > 1)
    {
      utl::getAttribute(elem,"order",phOrder);
      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        this->getIntegrand()->parse(child);
    }
    else if (!strcasecmp(elem->Value(),"elasticity") && !Dim::myProblem)
    {
      if (phOrder < 2)
      {
        std::string form("voigt");
        if (utl::getAttribute(elem,"formulation",form,true) && form != "voigt")
          Dim::myProblem = new FractureElasticity(Dim::dimension);
      }
      result = this->SIMElasticity<Dim>::parse(elem);
    }
    else
      result = this->SIMElasticity<Dim>::parse(elem);
    --ncall;
    return result;
  }

private:
  DynSIM dSim;    //!< Dynamic solution driver
  Matrix projSol; //!< Projected secondary solution fields
  Matrix eNorm;   //!< Element norm values
  Vector gNorm;   //!< Global norm values
  int    vtfStep; //!< VTF file step counter
  int    subIter; //!< Sub-iteration type flag
  int    phOrder; //!< Phase-field order for monolithic coupling

  std::string energFile; //!< File name for global energy output (monolithic)
};

#endif
