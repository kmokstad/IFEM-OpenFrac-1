// $Id$
//==============================================================================
//!
//! \file SIMFractureDynamics.h
//!
//! \date Jul 13 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Driver for fracture-dynamic problems.
//!
//==============================================================================

#ifndef _SIM_FRACTURE_DYNAMICS_H_
#define _SIM_FRACTURE_DYNAMICS_H_

#include "ProcessAdm.h"
#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#include "LRSpline/LRSplineSurface.h"
#endif
#include <fstream>
#include <numeric>


/*!
  \brief Driver class for fracture dynamics simulators.
  \details A fracture dynamics simulator is a coupling between
  a dynamic elasticity solver and a phase field solver.
*/

template<class SolidSolver, class PhaseSolver,
         template<class S1, class S2> class Coupling>
class SIMFracture : public Coupling<SolidSolver,PhaseSolver>
{
  //! Convenience type
  typedef Coupling<SolidSolver,PhaseSolver> CoupledSIM;

public:
  //! \brief The constructor initializes the references to the two solvers.
  SIMFracture(SolidSolver& s1, PhaseSolver& s2, const std::string& inputfile)
    : CoupledSIM(s1,s2), infile(inputfile), aMin{0.0} {}
  //! \brief Empty destructor.
  virtual ~SIMFracture() {}

  //! \brief Initializes and sets up field dependencies.
  virtual void setupDependencies()
  {
    this->S1.registerDependency(&this->S2,"phasefield",1);
    // The tensile energy is defined on integration points and not nodal points.
    // It is a global buffer array across all patches in the model.
    // Use an explicit call instead of normal couplings for this.
    this->S2.setTensileEnergy(this->S1.getTensileEnergy());
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp, bool firstS1 = true)
  {
    if (tp.step == 1 && this->S1.haveCrackPressure())
      // Start the initial step by solving the phase-field first
      if (!this->S2.solveStep(tp,false))
        return false;

    return this->CoupledSIM::solveStep(tp,firstS1);
  }

  //! \brief Saves the converged results to VTF-file of a given time step.
  //! \details It also writes global energy quantities to file for plotting.
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (!energFile.empty() && tp.step > 0 &&
        this->S1.getProcessAdm().getProcId() == 0)
    {
      std::ofstream os(energFile, tp.step == 1 ? std::ios::out : std::ios::app);

      Vector BF, RF;
      this->S1.getBoundaryForce(BF,this->S1.getSolutions(),tp);
      this->S1.getBoundaryReactions(RF);

      if (tp.step == 1)
      {
        size_t i;
        os <<"#t eps_e external_energy eps+ eps- eps_b |c|"
           <<" eps_d-eps_d(0) eps_d";
        for (i = 0; i < BF.size(); i++)
          os <<" load_"<< char('X'+i);
        for (i = 0; i < RF.size(); i++)
          os <<" react_"<< char('X'+i);
        os << std::endl;
      }

      os << std::setprecision(11) << std::setw(6) << std::scientific
         << tp.time.t;
      for (double n1 : this->S1.getGlobalNorms()) os <<" "<< n1;
      const Vector& n2 = this->S2.getGlobalNorms();
      os <<" "<< (n2.size() > 2 ? n2[1] : 0.0);
      os <<" "<< (n2.size() > 1 ? n2[n2.size()-2] : 0.0);
      os <<" "<< (n2.size() > 0 ? n2.back() : 0.0);
      for (double f : BF) os <<" "<< utl::trunc(f);
      for (double f : RF) os <<" "<< utl::trunc(f);
      os << std::endl;
    }

    return this->S2.saveStep(tp,nBlock) && this->S1.saveStep(tp,nBlock);
  }

  //! \brief Parses staggering parameters from an XML element.
  virtual void parseStaggering(const TiXmlElement*) {}

  //! \brief Assigns the file name for global energy output.
  void setEnergyFile(const char* fName)
  {
    if (fName)
    {
      energFile = fName;
      IFEM::cout <<"\tFile for global energy output: "<< energFile << std::endl;
    }
  }

  //! \brief Stores current solution state in an internal buffer.
  void saveState()
  {
    sols = this->S1.getSolutions();
    sols.push_back(this->S2.getSolution());
    hsol = this->S2.getHistoryField();
  }

  //! \brief Refines the mesh on the initial configuration.
  bool initialRefine(double beta, double min_frac, int nrefinements)
  {
    if (this->S2.getInitRefine() >= nrefinements)
      return true; // Grid is sufficiently refined during input parsing
    else if (this->S2.hasIC("phasefield"))
      return true; // No initial refinement when specified initial phase field

    TimeStep step0;
    int newElements = 1;
    for (step0.iter = 0; newElements > 0; step0.iter++)
      if (!this->S2.solveStep(step0))
        return false;
      else
        newElements = this->adaptMesh(beta,min_frac,nrefinements);

    return newElements == 0;
  }

  //! \brief Refines the mesh with transfer of solution onto the new mesh.
  int adaptMesh(double beta, double min_frac, int nrefinements)
  {
#ifdef HAS_LRSPLINE
    ASMu2D* pch = dynamic_cast<ASMu2D*>(this->S1.getPatch(1));
    if (!pch) return -999; // Logic error, should not happen

    if (aMin.front() <= 0.0) // maximum refinements per element
    {
      double redMax = pow(2.0,nrefinements);
      aMin.resize(this->S1.getNoPatches());
      for (int i = 0; i < this->S1.getNoPatches(); ++i) {
        pch = dynamic_cast<ASMu2D*>(this->S1.getPatch(i+1));
        aMin[i] = pch->getBasis()->getElement(0)->area()/(redMax*redMax);
      }
    }

    // Fetch element norms to use as refinement criteria
    Vector eNorm;
    double gNorm = this->S2.getNorm(eNorm,3);
    if (eNorm.empty())
    {
      std::cerr <<" *** SIMFractureDynamics:adaptMesh: Missing refinement"
                <<" indicators, expected as the 3rd element norm."<< std::endl;
      return -1;
    }

    // Sort element indices based on comparing values in eNorm
    IntVec idx(eNorm.size());
    std::iota(idx.begin(),idx.end(),0);
    std::sort(idx.begin(),idx.end(),
              [&eNorm](size_t i1, size_t i2) { return eNorm[i1] < eNorm[i2]; });

    double eMin = min_frac < 0.0 ? -min_frac*gNorm/sqrt(idx.size()) : min_frac;
    size_t eMax = beta < 0.0 ? idx.size() : idx.size()*beta/100.0;
    IFEM::cout <<"\n  Lowest element: "<< std::setw(8) << idx.front()
               <<"    |c| = "<< eNorm[idx.front()]
               <<"\n  Highest element:"<< std::setw(8) << idx.back()
               <<"    |c| = "<< eNorm[idx.back()]
               <<"\n  Minimum |c|-value for refinement: "<< eMin
               <<"\n  Minimum element area:";
    for (double amin : aMin)
      IFEM::cout << " " << amin;
    IFEM::cout << std::endl;

    IntVec elements; // Find the elements to refine
    elements.reserve(eMax);
    for (int eid : idx)
      if (eNorm[eid] > eMin || elements.size() >= eMax)
        break;
      else {
        size_t p = 0;
        for (const ASMbase* patch : this->S1.getFEModel()) {
          int el;
          ++p;
          if ((el = patch->getElmIndex(eid+1))) {
            const ASMu2D* pch = dynamic_cast<const ASMu2D*>(patch);
            if (pch->getBasis()->getElement(el-1)->area() > aMin[p-1]+1.0e-12)
              elements.push_back(eid);
            break;
          }
        }
      }

    if (elements.empty())
      return 0;

    IFEM::cout <<"  Elements to refine: "<< elements.size()
               <<" (|c| = ["<< eNorm[elements.front()]
               <<","<< eNorm[elements.back()] <<"])\n"<< std::endl;

    std::vector<LR::LRSplineSurface*> oldBasis;
    if (!hsol.empty())
      for (ASMbase* patch : this->S1.getFEModel()) {
        ASMu2D* pch = dynamic_cast<ASMu2D*>(patch);
        oldBasis.push_back(pch->getBasis()->copy());
      }

    // Do the mesh refinement
    LR::RefineData prm;
    prm.options = { 10, 1, 2, 0, this->S1.getNoPatches() > 1 ? -1 : 1 };

    // Translate from element IDs to function IDs
    for (int elem : elements) {
      for (const ASMbase* patch : this->S1.getFEModel()) {
        int el;
        if ((el = patch->getElmIndex(elem+1))) {
          const ASMu2D* pch = dynamic_cast<const ASMu2D*>(patch);
          for (const LR::Basisfunction* b : pch->getBasis()->getElement(el-1)->support())
            prm.elements.push_back(pch->getNodeID(b->getId()+1)-1);
          break;
        }
      }
    }

    size_t solsize = sols.size();
    // sanity for some case. cannot recall
    size_t sollen = sols.empty() ? 0 : sols.front().size();
    if (!this->S1.refine(prm,sols) || !this->S2.refine(prm))
      return -2;

    // Re-initialize the simulators for the new mesh
    this->S1.clearProperties();
    this->S2.clearProperties();
    if (!this->S1.read(infile.c_str()) || !this->S2.read(infile.c_str()))
      return -3;

    if (!this->preprocess())
      return -4;

    if (!this->init(TimeStep()))
      return -5;

    if (!this->S1.initSystem(this->S1.opt.solver,1,1,0,true) ||
        !this->S2.initSystem(this->S2.opt.solver))
      return -6;

    // Transfer solution variables onto the new mesh
    if (!sols.empty())
    {
      IFEM::cout <<"\nTransferring "<< solsize-1 <<"x"<< sollen
                 <<" solution variables to new mesh for "<< this->S1.getName();
      Vectors soli(solsize-1, Vector(this->S1.getNoDOFs()));
      for (size_t i = 0; i < solsize-1; ++i) {
        for (int p = 0; p < this->S1.getNoPatches(); ++p)
          this->S1.injectPatchSolution(soli[i], sols[p*solsize+i], p);
      }
      this->S1.setSolutions(soli);
      IFEM::cout <<"\nTransferring "<< sollen
                 <<" solution variables to new mesh for "<< this->S2.getName();
      Vector solc(this->S2.getNoDOFs());
      for (int p = 0; p < this->S1.getNoPatches(); ++p)
        this->S2.injectPatchSolution(solc, sols[p*solsize+solsize-1], p);
      this->S2.setSolution(solc);
    }
    if (!hsol.empty())
    {
      IFEM::cout <<"\nTransferring "<< hsol.size()
                 <<" history variables to new mesh for "<< this->S2.getName()
                 << std::endl;
      this->S2.transferHistory2D(hsol,oldBasis);
      for (LR::LRSplineSurface* srf : oldBasis)
        delete srf;
    }

    return elements.size();
#else
    std::cerr <<" *** SIMFractureDynamics:adaptMesh: No LR-spline support.\n";
    return -1;
#endif
  }

private:
  std::string energFile; //!< File name for global energy output
  std::string infile;    //!< Input file parsed

  RealArray aMin; //!< Minimum element area
  Vectors   sols; //!< Solution state to transfer onto refined mesh
  RealArray hsol; //!< History field to transfer onto refined mesh
};

#endif
