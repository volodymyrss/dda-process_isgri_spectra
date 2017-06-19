import ddosa 
from astropy.io import fits 
from bcolors import render
import subprocess
import os
import dataanalysis as da
from numpy import *

try:
    import crab
except:
    pass

def get_open_fds():
    '''
    return the number of open file descriptors for current process

    .. warning: will only work on UNIX-like os-es.
    '''

    pid = os.getpid()
    procs = subprocess.check_output( 
        #[ "lsof", '-w', '-Ff', "-p", str( pid ) ] )
        [ "lsof", '-w', '-Fn', "-p", str( pid ) ] )

    #files=filter( 
    #        lambda s: s and s[ 0 ] == 'f' and s[1: ].isdigit(),
    files=        procs.split( '\n' ) 

    #print files

    nprocs = len( 
            files
        )
    return nprocs

class ProcessSpectra(ddosa.DataAnalysis):
    input_spectra=ddosa.ii_spectra_extract
    #response="/resources/rmf_256bands.fits"
    #arf="/Integral/data/resources/arfs/arf_62_1528.fits"
    arf=None

    def main(self):
        f=fits.open(self.input_spectra.spectrum.path)


        for hdu in f[2:]:
            sname=hdu.header['NAME']
            fn="isgri_spectrum_%s.fits"%sname.replace(" ","_")
       #     hdu.header['RESPFILE']=self.response

            if self.arf is not None:
                hdu.header['ANCRFILE']=self.arf
            print "source:",sname,"to",fn
            print hdu.header['RESPFILE']
            print hdu.header['ANCRFILE']
            hdu.writeto(fn,clobber=True)
            setattr(self,fn,da.DataFile(fn))




class ScWSpectraList(ddosa.DataAnalysis):
    input_scwlist=ddosa.RevScWList
    copy_cached_input=False
    input_specsummary=ddosa.SpectraProcessingSummary

    allow_alias=True

    version="allthem"
    
    maxspec=None

    def main(self):
        self.spectra=[ddosa.ii_spectra_extract(assume=scw) for scw in self.input_scwlist.scwlistdata]
        
        if self.maxspec is not None: self.spectra=self.spectra[:self.maxspec]

class SpectrumFromFile(da.DataAnalysis):
    input_filename=None
        
    def main(self):
        self.spectrum=da.DataFile(self.input_filename.str())

class FileSpectraList(ddosa.DataAnalysis):
    input_file=None
    copy_cached_input=False

    allow_alias=True

    version="allthem"
    
    maxspec=None

    def main(self):
        self.spectra=[SpectrumFromFile(input_filename=fn) for fn in open(self.input_file.filename)]

        print self.spectra
        
        if self.maxspec is not None: self.spectra=self.spectra[:self.maxspec]

class SpectrumEfficiencyCorrection(ddosa.DataAnalysis):
    input_bins=ddosa.SpectraBins

    enable=True

    def main(self):
        print "will do"

    def get_version(self):
        v=self.get_signature()+"."+self.version
        if not self.enable:
            v+=".disabled"
        else:
            v+=".enabled"
        return v

    def correct(self,r,e):
        if not self.enable: return r,e

        e1,e2=map(array,zip(*self.input_bins.bins))
        ec=(e1+e2)/2.

        a = 13.0435898895737
        b = -6.81616014518407
        c = 0.87794025777641
        f=lambda x:exp(a+log(x)*b+log(x)**2*c)

        cr=r[:]
        m=(ec>30) & (ec<80)
        cr[m]/=f(ec)[m]

        return cr,e

class ISGRISpectraSum(ddosa.DataAnalysis):
    input_spectralist=ScWSpectraList

    input_response=ddosa.SpectraBins

    input_efficiency=SpectrumEfficiencyCorrection

    copy_cached_input=False


    cached=True

    version="v5.1"

    sources=['Crab']

    extract_all=False
    save_lc=True

    def main(self):
        spectra={}

        choice=self.input_spectralist.spectra

        allsource_summary=[]


        sig=lambda x,y:(((x/y)[~isnan(y) & (y!=0) & ~isinf(y) & ~isinf(x) & ~isnan(x)])**2).sum()**0.5

        import time
        t0=time.time()
        i_spec=1

        for spectrum in choice:
            if hasattr(spectrum,'empty_results'):
                print "skipping",spectrum
                continue

            fn=spectrum.spectrum.get_path()
            print "%i/%i"%(i_spec,len(choice))
            tc=time.time()
            print "seconds per spectrum:",(tc-t0)/i_spec,"will be ready in %.5lg seconds"%((len(choice)-i_spec)*(tc-t0)/i_spec)
            i_spec+=1
            print "spectrum from",fn

            f=fits.open(fn)


            t1,t2=f[1].header['TSTART'],f[1].header['TSTOP']
            print t1,t2

            for e in f:
                try:
                    if e.header['EXTNAME']!="ISGR-EVTS-SPE": continue
                except:
                    continue

                try:
                    name=e.header['NAME']
                except:
                    name="Unnamed"
                allsource_summary.append([name,t1,t2,copy(e.data['RATE']),copy(e.data['STAT_ERR'])])
                if (name in self.sources) or (self.extract_all):
                    rate=e.data['RATE']
                    err=e.data['STAT_ERR']
                    exposure=e.header['EXPOSURE']
                    if name not in spectra:
                        spectra[name]=[rate,err**2,exposure,e]
                        preserve_file=True
                    else:
                        err[isnan(err) | (err==0)]=inf
                        spectra[name][1][isnan(spectra[name][1]) | (spectra[name][1]==0)]=inf
                        spectra[name][0]=(spectra[name][0]/spectra[name][1]+rate/err**2)/(1/spectra[name][1]+1/err**2)
                        #spectra[name][0]=(spectra[name][0]/spectra[name][1]+rate/err**2)/(1/spectra[name][1]+1/err**2)
                        spectra[name][1]=1/(1/spectra[name][1]+1/err**2)
                        spectra[name][2]+=exposure
            
                    print render("{BLUE}%.20s{/}"%name),"%.4lg sigma in %.5lg ks"%(sig(rate,err),exposure/1e3),"total %.4lg in %.5lg ks"%(sig(spectra[name][0],spectra[name][1]), spectra[name][2]/1e3)

            #if not preserve_file:
            #    print "closing file"
            f.close()

            try:
                print get_open_fds()
            except Exception as e:
                print "unable to check open fds"

        eb1,eb2=map(array,zip(*self.input_response.bins))

        source_results=[]

        

        for name,spectrum in spectra.items():
            source_short_name=name.strip().replace(" ","_")

            rmf_fn="rmf_sum_%s.fits"%source_short_name
            fits.open(self.input_response.binrmf).writeto(rmf_fn,clobber=True)
            #fits.open(spectrum[3].header['RESPFILE']).write(rmf_fn,clobber=True)
            
            if spectrum[3].header['ANCRFILE'] == 'NONE':
                arf_fn=None
            else:
                arf_fn="arf_sum_%s.fits"%source_short_name

                try:
                    fits.open(spectrum[3].header['ANCRFILE']).writeto(arf_fn,clobber=True)
                except IOError:
                    arf_fn=None

            spectrum[3].data['RATE'][:],spectrum[3].data['STAT_ERR'][:]=self.input_efficiency.correct(spectrum[0][:],(spectrum[1]**0.5)[:])
            spectrum[3].header['EXPOSURE']=spectrum[2]
            #spectrum[3].header['RESPFILE']=self.input_response.binrmf

            spectrum[3].header['RESPFILE']=rmf_fn
            spectrum[3].header['ANCRFILE']=arf_fn

            fn="isgri_sum_%s.fits"%source_short_name
            spectrum[3].writeto(fn,clobber=True)


            select_range=lambda x,a,b:((eb1>a) & (eb2<b) & ~isnan(x) & ~isinf(x))

            print render("{RED}total{/} for {BLUE}%s{/}"%name)

            all_spectra=[l for l in allsource_summary if l[0]==name]
            #print all_spectra

            for erange in [(25,80),(80,300),(10,800)]:
                m=select_range(spectrum[3].data['RATE'],erange[0],erange[1])
                r=spectrum[3].data['RATE'][m]
                e=spectrum[3].data['STAT_ERR'][m]

                lc=[(l[1],l[2],sum(l[3][m]),sum(l[4][m]**2)**0.5) for l in all_spectra]

                lc_t1,lc_t2,lc_f,lc_fe=map(array,zip(*lc))

                if self.save_lc:
                    savetxt("%s_%.5lg_%.5lg.txt"%(name.replace(" ","_"),erange[0],erange[1]),column_stack((lc_t1,lc_t2,lc_f,lc_fe)))

                varamp=std(lc_f)
                #varfrac=std(lc_f)/average(lc_fe**2)**0.5
                
                try:
                    mCrab=crab.Crab().counts_in(erange[0],erange[1])/1000.
                except:
                    mCrab=1000000
                source_stats=[erange[0],erange[1],r.sum(),(e**2).sum()**0.5,r.sum()/mCrab,sig(r,e),varamp]
                print render("{RED}total{/} %.5lg-%.5lg keV %.5lg+-%.5lg cts %.5lg mCrab %.5lg sigma, variability %.5lg"%tuple(source_stats))

                source_results.append([name]+source_stats)
                

            setattr(self,fn.replace(".fits",""),da.DataFile(fn))
            
            if rmf_fn is not None:
                setattr(self,rmf_fn.replace(".fits",""),da.DataFile(rmf_fn))

            if arf_fn is not None:
                setattr(self,arf_fn.replace(".fits",""),da.DataFile(arf_fn))

        srf=open("source_summary.txt","w")
        for sr in source_results:
            srf.write(sr[0].replace(" ","_")+" "+" ".join(["%.5lg"%s for s in sr[1:]])+"\n")
            
 #       for l in allsource_summary:
#            print l[0] #,l[1].shape
