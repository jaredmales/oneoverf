#include <vector>
#include <cmath>

#include <mx/app/application.hpp>

#include <mx/improc/eigenImage.hpp>
#include <mx/improc/eigenCube.hpp>
#include <mx/ioutils/fits/fitsFile.hpp>
#include <mx/sigproc/psdUtils.hpp>
#include <mx/sigproc/psdFilter.hpp>

#include <mx/math/randomT.hpp>
#include <mx/math/vectorUtils.hpp>

#include <mx/ipc/ompLoopWatcher.hpp>

namespace mx
{
namespace sigproc
{


}
}

/// A program to calculate correlated noise 3D cubes
/**
  * 
  */
template<typename _realT>
class oneoverf : public mx::app::application
{
public:
   typedef _realT realT;


protected:
   
   realT m_space_alpha {2};
   realT m_space_L0 {0};
   realT m_space_l0 {0};
   int   m_space_dim {128};
   int   m_space_oversamp {2};

   realT m_time_alpha {2};
   realT m_time_T0 {0};
   realT m_time_t0 {0};
   int   m_time_dim {64};
   int   m_time_oversamp {2};

   realT m_rms {1};
   bool m_force_rms {true};

   std::string m_outputName {"filtered.fits"};

public:

   oneoverf();

   ~oneoverf();

   /// Setup the configuration system
   /**
     */
   virtual void setupConfig();

   /// Load the configuration
   /**
     */
   virtual void loadConfig();

   /// Do the calculations.
   /**
     */
   virtual int execute();

};

template<typename realT>
oneoverf<realT>::oneoverf()
{
}

template<typename realT>
oneoverf<realT>::~oneoverf()
{
}

template<typename realT>
void oneoverf<realT>::setupConfig()
{
   config.add("space.alpha",    "", "space.alpha",    mx::app::argType::Required, "space", "alpha",    false, "real", "space (2D) power law exopnent.  Default 2 (1/f^2).");
   config.add("space.L0",       "", "space.L0",       mx::app::argType::Required, "space", "L0",       false, "real", "space (2D) outer scale.  Default 0 (no outer scale)");
   config.add("space.l0",       "", "space.l0",       mx::app::argType::Required, "space", "l0",       false, "real", "space (2D) inner scale.  Default 0 (no inner scale)");
   config.add("space.dim",      "", "space.dim",      mx::app::argType::Required, "space", "dim",      false, "int",  "space (2D) dimension.  Images will be dim x dim.  Default 128.");
   config.add("space.oversamp", "", "space.oversamp", mx::app::argType::Required, "space", "oversamp", false, "int",  "space (2D) oversampling.  Default 2");

   config.add("time.alpha",    "", "time.alpha",    mx::app::argType::Required, "time", "alpha",    false, "real", "time (1D) power law exopnent.  Default 2 (1/f^2).");
   config.add("time.T0",       "", "time.T0",       mx::app::argType::Required, "time", "T0",       false, "real", "time (1D) outer scale.  Default 0 (no outer scale)");
   config.add("time.t0",       "", "time.t0",       mx::app::argType::Required, "time", "t0",       false, "real", "time (1D) inner scale.  Default 0 (no inner scale)");
   config.add("time.dim",      "", "time.dim",      mx::app::argType::Required, "time", "dim",      false, "int",  "time (1D) dimension.  Cube will contain dim images.  Default 1024.");
   config.add("time.oversamp", "", "time.oversamp", mx::app::argType::Required, "time", "oversamp", false, "int",  "time (1D) oversampling.  Default 2");

   config.add("rms.target",      "", "rms.target", mx::app::argType::Required, "rms", "target", false, "real", "Target rms.  Default 1");
   config.add("rms.force",       "", "rms.force",  mx::app::argType::Required, "rms", "force",  false, "bool", "If true, mean is set to 0 and output scaled to have rms.target.  If false, then rms will only be the ensemble value.  Default true");

   config.add("output.name",  "", "output.name",  mx::app::argType::Required, "output", "name",  false, "string", "Output file name. Default filtered.fits");
}

template<typename realT>
void oneoverf<realT>::loadConfig()
{ 
   config(m_space_alpha, "space.alpha");
   config(m_space_L0, "space.L0");
   config(m_space_l0, "space.l0");
   config(m_space_dim, "space.dim");
   config(m_space_oversamp, "space.oversamp");

   config(m_time_alpha, "time.alpha");
   config(m_time_T0, "time.T0");
   config(m_time_t0, "time.t0");
   config(m_time_dim, "time.dim");
   config(m_time_oversamp, "time.oversamp");


   config(m_rms, "rms.target");
   config(m_force_rms, "rms.force");

   config(m_outputName, "output.name");

}

template<typename realT>
int oneoverf<realT>::execute()
{
   if(m_space_oversamp < 1)
   {
      std::cerr << "Setting space.oversamp to 1\n";
      m_space_oversamp = 1;
   }

   if(m_time_oversamp < 1)
   {
      std::cerr << "Setting time.oversamp to 1\n";
      m_time_oversamp = 1;
   }

   //---- make the temporal PSD ----
   std::vector<realT> f, tpsd;
   f.resize(m_time_dim * m_time_oversamp);
   mx::sigproc::frequencyGrid(f, 1.0);

   tpsd.resize(f.size());

   mx::sigproc::vonKarmanPSD(tpsd, f, (realT) 1.0, m_time_T0, m_time_t0, m_time_alpha);
   if( !std::isnormal(tpsd[0]) || std::isinf(tpsd[0])) tpsd[0] = 0; //Make sure we don't propagate a nan/inf

   mx::sigproc::normPSD(tpsd, f, 1.0); //normalize it

   // std::cerr << fmin << "\n";
   // for(size_t n=0; n< tpsd.size(); ++n)
   // {
   //    std::cout << f[n] << " " << tpsd[n] << "\n";
   // }
   // std::cerr << mx::math::vectorSum(tpsd)*f[1] << "\n";

   // exit(0);

   //---- allocate the oversized phase screen cube ----
   std::cerr << "allocating...\n";
   mx::improc::eigenCube<realT> tims;
   tims.resize( m_space_dim*m_space_oversamp, m_space_dim*m_space_oversamp, m_time_dim);

   //---- fill the temporally filtered cube ----
   mx::ipc::ompLoopWatcher<std::ostream> tflw(tims.cols(), std::cerr);

   std::cerr << "temporal filtering...\n";
   #pragma omp parallel
   {
      mx::math::normDistT<realT> normVar; //produces standard normal variates
      
      //The temporary time series for filtering
      std::vector<realT> tts;
      tts.resize(m_time_dim * m_time_oversamp);

      //The 1D PSD filter
      mx::sigproc::psdFilter<realT,1> tfilt;
      tfilt.psd(tpsd,f[1]-f[0]);
   
      #pragma omp for
      for(int c=0; c < tims.cols(); ++c)
      {
         for(int r=0; r < tims.rows(); ++r)
         {
            //Fill in oversampled white noise time series
            for(size_t n =0; n < tts.size(); ++n)
            {
               tts[n] = normVar;
            }

            //filter it
            tfilt(tts);

            //Now extract the sample 
            for(size_t n=0; n < m_time_dim; ++n)
            {
               tims.image(n)(r,c) = tts[n];
            }
         }
         tflw.incrementAndOutputStatus();

      }//for c
   }//pragma omp parallel

   //---- make the spatial PSD ----
   mx::improc::eigenImage<realT> k, spsd;
   k.resize( m_space_dim*m_space_oversamp, m_space_dim*m_space_oversamp);
   mx::sigproc::frequencyGrid(k, 1.0);
   
   spsd.resize(k.rows(), k.cols());
   mx::sigproc::vonKarman_psd(spsd, k, m_space_alpha, m_space_L0, m_space_l0);

   mx::sigproc::normPSD(spsd, k, pow(m_rms,2)); //normalize

   //---- allocate the final phase screen cube ----
   mx::improc::eigenCube<realT> ims;
   ims.resize( m_space_dim, m_space_dim, m_time_dim);

   //---- spatially filter the cube ----
   mx::ipc::ompLoopWatcher<std::ostream> sflw(tims.planes(), std::cerr);

   std::cerr << "spatial filtering...                                            \n";
   #pragma omp parallel
   {
      //the 2D filter
      mx::sigproc::psdFilter<realT,2> sfilt;
      sfilt.psd(spsd,k(1,0)-k(0,0), k(1,0)-k(0,0));
   
      mx::improc::eigenImage<realT> im;

      #pragma omp for
      for(int p=0; p<tims.planes(); ++p)
      {
         sfilt.filter(tims.image(p));

         for(int c=0; c < ims.rows(); ++c)
         {
            for(int r=0; r < ims.cols(); ++r)
            {
               ims.image(p)(r,c) = tims.image(p)(r,c);
            }
         }

         sflw.incrementAndOutputStatus();
      }
   }

   //tims.clear(); //let the memory go!
   realT mn = mx::math::vectorMean(ims.data(), ims.planes()*ims.rows()*ims.cols());
   realT var, rms;
   if(m_force_rms)
   {
      var = mx::math::vectorVariance(ims.data(), ims.planes()*ims.rows()*ims.cols(), mn);

      rms = sqrt(var);

      for(size_t n=0 ; n < ims.planes()*ims.rows()*ims.cols(); ++n) ims.data()[n] = (ims.data()[n] - mn)*(m_rms/rms);
   }

   var = mx::math::vectorVariance(ims.data(), ims.planes()*ims.rows()*ims.cols(), (realT) 0.0);
   rms = sqrt(var);

   mn = mx::math::vectorMean(ims.data(), ims.planes()*ims.rows()*ims.cols());

   std::cerr << "Mean = " << mn << "                                           \n";
   std::cerr << "rms  = " << rms << "\n";

   mx::fits::fitsHeader fh;

   fh.append("TMALPHA", m_space_alpha, "Time PSD exponent alpha");
   fh.append("TMOUTSCL", m_space_L0, "Time PSD outer scale T0");
   fh.append("TMINNSCL", m_space_l0, "Time PSD inner scale t0");
   fh.append("TMOVRSMP", m_space_oversamp, "Time oversampling");

   fh.append("SPALPHA", m_space_alpha, "Space PSD exponent alpha");
   fh.append("SPOUTSCL", m_space_L0, "Space PSD outer scale L0");
   fh.append("SPINNSCL", m_space_l0, "Space PSD inner scale l0");
   fh.append("SPOVRSMP", m_space_oversamp, "Space oversampling");

   fh.append("RMSTGT", m_rms, "Target rms");
   fh.append("RMSFORCE", m_force_rms, "RMS force flag");
   fh.append("MEAN", mn, "Measured mean");
   fh.append("RMS", rms, "Measured RMS");

   mx::fits::fitsFile<realT> ff;

   ff.write(m_outputName, ims, fh);

   std::cout << "written to: " << m_outputName << "\n";
   
   return 0;
}


int main( int argc,
          char ** argv
        )
{
   oneoverf<float> oof;

   return oof.main(argc,argv);
}

