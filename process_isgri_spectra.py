from __future__ import print_function

import ddosa 
from astropy.io import fits 
from pscolors import render
import datetime
import subprocess
import shutil
import os

try:
    from dataanalysis import core as da
except ImportError:
    import dataanalysis as da

from numpy import *
from collections import defaultdict

import useresponse

#class MultiEpochRMFNotSupported(da.AnalysisException):
class MultiEpochRMFNotSupported(Exception):
    pass

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

    #print(files)

    nprocs = len( 
            files
        )
    return nprocs

class ProcessSpectra(ddosa.DataAnalysis):
    input_spectra=ddosa.ii_spectra_extract
    input_arf=useresponse.FindResponse
    input_response=ddosa.ISGRIResponse


    def main(self):
        f=fits.open(self.input_spectra.spectrum.path)


        for hdu in f[2:]:
            sname=hdu.header['NAME']
            fn="isgri_spectrum_%s.fits"%sname.replace(" ","_")

            hdu.header['ANCRFILE']=self.input_arf.arf_path
            hdu.header['RESPFILE']=self.input_response.path

            print("source:",sname,"to",fn)
            print(hdu.header['RESPFILE'])
            print(hdu.header['ANCRFILE'])
            hdu.writeto(fn,overwrite=True)
            setattr(self,fn,da.DataFile(fn))


class ISGRISpectrumPack(ddosa.DataAnalysis):
    input_spectra=ddosa.ii_spectra_extract
    input_arf=useresponse.FindResponse
    input_response=ddosa.ISGRIResponse

    cached=True

    def main(self):
        if hasattr(self.input_spectra,'empty_results'):
            print("skipping")
            self.empty_results = True
            return

    
        self.spectra_spectrum = ddosa.localized_DataFile(self.input_spectra.spectrum)
        self.arf = ddosa.localized_DataFile(self.input_arf.arf_path)
        self.rmf = ddosa.localized_DataFile(self.input_response.path)


class ScWSpectraList(ddosa.DataAnalysis):
    input_scwlist=ddosa.RevScWList
    copy_cached_input=False
    input_specsummary=ddosa.SpectraProcessingSummary

    allow_alias=True

    version="allthem"
    
    maxspec=None

    def main(self):
        self.spectra=[ISGRISpectrumPack(assume=scw) for scw in self.input_scwlist.scwlistdata]
        #self.spectra=[[ddosa.ii_spectra_extract(assume=scw),useresponse.FindResponse(assume=scw),ddosa.ISGRIResponse(assume=scw)] for scw in self.input_scwlist.scwlistdata]

        if len(self.spectra)==0:
            raise ddosa.EmptyScWList()
        
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
        self.spectra=[SpectrumFromFile(input_filename=fn) for fn in open(self.input_file.filename)] # incompat!

        print(self.spectra)
        
        if self.maxspec is not None: self.spectra=self.spectra[:self.maxspec]

class SpectrumEfficiencyCorrection(ddosa.DataAnalysis):
    input_bins=ddosa.SpectraBins

    enable=False

    def main(self):
        print("will do")

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


    #input_efficiency=SpectrumEfficiencyCorrection

    copy_cached_input=False

    spectra=None

    cached=True

    version="v5.4.2.4"

    sources=['Crab']

    extract_all=True
    save_lc=True

    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.extract_all:
            v+=".extractall"
        else:
            v+".extract_"+("_".join(self.sources))
        return v

    def main(self):
        spectra={}

        choice=self.input_spectralist.spectra

        allsource_summary=[]


        sig=lambda x,y:(((x/y)[~isnan(y) & (y!=0) & ~isinf(y) & ~isinf(x) & ~isnan(x)])**2).sum()**0.5

        import time
        t0=time.time()
        i_spec=1        

        used_spectra = []

        #for spectrum,arf,rmf in choice:
        for pack in choice:
            if hasattr(pack,'empty_results'):
                print("skipping", pack)
                continue
            
            if not hasattr(pack,'spectra_spectrum'):
                print("skipping", pack)
                continue

            fn=pack.spectra_spectrum.get_path()
            print("%i/%i"%(i_spec,len(choice)))
            tc=time.time()
            print("seconds per spectrum:",(tc-t0)/i_spec,"will be ready in %.5lg seconds"%((len(choice)-i_spec)*(tc-t0)/i_spec))
            i_spec+=1
            print("spectrum from",fn)

            f=fits.open(fn).copy()

            if fn in used_spectra:
                raise Exception("identical spectra in different scw! impossible!")
            else:
                used_spectra.append(fn)

            try:
                t1,t2=f[1].header['TSTART'],f[1].header['TSTOP']
            except IndexError as e:                
                print("skipping for no start/stop", pack)
                continue

            print(t1,t2)

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
                    exp_src=e.header['EXP_SRC']
                    ontime=e.header['ONTIME']
                    telapse=e.header['TELAPSE']
                    tstart=e.header['TSTART']
                    tstop=e.header['TSTOP']
                    revol=e.header['REVOL']
                    if name not in spectra:
                        spectra[name]=[rate,err**2,exposure,e,defaultdict(int),defaultdict(int),ontime,telapse,tstart,tstop,[revol],exp_src]
                        preserve_file=True
                    else:
                        err[isnan(err) | (err==0)]=inf
                        spectra[name][1][isnan(spectra[name][1]) | (spectra[name][1]==0)]=inf
                        spectra[name][0]=(spectra[name][0]/spectra[name][1]+rate/err**2)/(1/spectra[name][1]+1/err**2)
                        #spectra[name][0]=(spectra[name][0]/spectra[name][1]+rate/err**2)/(1/spectra[name][1]+1/err**2)
                        spectra[name][1]=1/(1/spectra[name][1]+1/err**2)
                        spectra[name][2]+=exposure
                        spectra[name][6]+=ontime
                        spectra[name][7]+=telapse
                        spectra[name][8]=min(spectra[name][8],tstart)
                        spectra[name][9]=max(spectra[name][9],tstop)
                        if revol not in spectra[name][10]:
                            spectra[name][10].append(revol)
                        spectra[name][11]+=exp_src

                    arf_path=pack.arf.get_path()
                    rmf_path=pack.rmf.get_path()

                    spectra[name][4][arf_path]+=exposure
                    spectra[name][5][rmf_path]+=exposure
            
                    print(render("{BLUE}%.20s{/}"%name),"%.4lg sigma in %.5lg ks"%(sig(rate,err),exposure/1e3),"total %.4lg in %.5lg ks"%(sig(spectra[name][0],spectra[name][1]), spectra[name][2]/1e3))


            #if not preserve_file:
            #    print("closing file")
            f.close()

            try:
                print(get_open_fds())
            except Exception as e:
                print("unable to check open fds")

        eb1,eb2=map(array,zip(*self.input_response.bins))

        self.spectra=spectra

        source_results=[]
        self.extracted_sources=[]

        for name,spectrum in spectra.items():
            source_short_name=name.strip().replace(" ","_").replace("/","")

           # if len(spectrum[5].keys())>1 or len(spectrum[4].keys())>1:
           #     raise MultiEpochNotSupported()

            # sum arf

            arfs=sorted(spectrum[4].items())
            print("got",len(arfs),"arfs:")
            for arf_fn,arf_exposure in arfs:
                print(arf_fn,arf_exposure)

            print("summing arfs")
            arf_first=fits.open(arfs[0][0])
            arf_first[1].data['SPECRESP']*=arfs[0][1]
            total_exposure=arfs[0][1]
            print("arf:",arfs[0][0],arf_first[1].data['SPECRESP'].max(),"cm**2 * s")

            for arf_fn,arf_exposure in arfs[1:]:
                arf_f=fits.open(arf_fn)
                arf_first[1].data['SPECRESP']+=arf_f[1].data['SPECRESP']*arf_exposure
                total_exposure+=arf_exposure
                print("arf:",arf_fn,arf_f[1].data['SPECRESP'].max()*arf_exposure,"cm**2 * s")
            arf_first[1].data['SPECRESP']/=total_exposure
                

            rmf_checksums = []
            for rmf_fn, rmf_exposure in spectrum[5].items():
                cs = fits.open(rmf_fn)[1].header['CHECKSUM']
                rmf_checksums.append(cs)

            if len(list(set(rmf_checksums))) > 1:
                raise MultiEpochRMFNotSupported(f"found many epochs: {list(spectrum[5].keys())}")

            arf_fn="arf_sum_%s.fits"%source_short_name
            #fits.open(spectrum[4].keys()[0]).writeto(arf_fn,overwrite=True)
            print("total arf:",arf_fn,arf_first[1].data['SPECRESP'].max()*total_exposure,"cm**2 * s")
            arf_first.writeto(arf_fn,overwrite=True)
            

            rmf_fn="rmf_sum_%s.fits"%source_short_name
            print("individual rmf",list(spectrum[5].keys())[0])
            fits.open(list(spectrum[5].keys())[0]).writeto(rmf_fn,overwrite=True)
            

            #spectrum[3].data['RATE'][:],spectrum[3].data['STAT_ERR'][:]=self.input_efficiency.correct(spectrum[0][:],(spectrum[1]**0.5)[:])
            spectrum[3].data['RATE'][:],spectrum[3].data['STAT_ERR'][:]=spectrum[0][:],(spectrum[1]**0.5)[:]
            spectrum[3].header['EXPOSURE']=spectrum[2]
            spectrum[3].header['ONTIME']=spectrum[6]
            spectrum[3].header['TELAPSE']=spectrum[7]
            spectrum[3].header['EXP_SRC']=spectrum[11]
            #spectrum[3].header['RESPFILE']=self.input_response.binrmf

            spectrum[3].header['RESPFILE']=rmf_fn
            spectrum[3].header['ANCRFILE']=arf_fn

            spectrum[3].header['CREATOR']=self.get_version()
            spectrum[3].header['DATE']=datetime.datetime.now().isoformat()

            if len(spectrum[10])==1:
                spectrum[3].header['REVOL']=spectrum[10][0]
            else:
                del spectrum[3].header['REVOL']
                if len(spectrum[10])<4:
                    spectrum[3].header['REVOLS']=",".join("%.4i"%r for r in spectrum[10])
                else:
                    spectrum[3].header['REVOLS']="%.4i..%.4i"%(sorted(spectrum[10])[0],sorted(spectrum[10])[-1])

            del spectrum[3].header['SWID']
            del spectrum[3].header['SWBOUND']
            del spectrum[3].header['OBTSTART']
            del spectrum[3].header['OBTEND']
            spectrum[3].header['TSTART']=spectrum[8]
            spectrum[3].header['TSTOP']=spectrum[9]
            del spectrum[3].header['TFIRST']
            del spectrum[3].header['TLAST']
            del spectrum[3].header['Y_FIN']
            del spectrum[3].header['Z_FIN']
            spectrum[3].header['STAMP']=None


            fn="isgri_sum_%s.fits"%source_short_name
            spectrum[3].writeto(fn, overwrite=True, checksum=True)
            print("writing",fn)


            select_range=lambda x,a,b:((eb1>a) & (eb2<b) & ~isnan(x) & ~isinf(x))

            print(render("{RED}total{/} for {BLUE}%s{/}"%name))

            all_spectra=[l for l in allsource_summary if l[0]==name]
            #print(all_spectra)

            for erange in [(25,80),(80,300),(10,800)]:
                m=select_range(spectrum[3].data['RATE'],erange[0],erange[1])
                r=spectrum[3].data['RATE'][m]
                e=spectrum[3].data['STAT_ERR'][m]

                lc=[(l[1],l[2],sum(l[3][m]),sum(l[4][m]**2)**0.5) for l in all_spectra]

                lc_t1,lc_t2,lc_f,lc_fe=map(array,zip(*lc))

                if self.save_lc:
                    lc_fn="%s_%.5lg_%.5lg.txt"%(source_short_name,erange[0],erange[1])
                    savetxt(lc_fn,column_stack((lc_t1,lc_t2,lc_f,lc_fe)))
                    setattr(self,lc_fn.replace(".txt",""),da.DataFile(lc_fn))

                varamp=std(lc_f)
                #varfrac=std(lc_f)/average(lc_fe**2)**0.5
                
                try:
                    mCrab=crab.Crab().counts_in(erange[0],erange[1])/1000.
                except:
                    mCrab=1000000
                source_stats=[erange[0],erange[1],r.sum(),(e**2).sum()**0.5,r.sum()/mCrab,sig(r,e),varamp]
                print(render("{RED}total{/} %.5lg-%.5lg keV %.5lg+-%.5lg cts %.5lg mCrab %.5lg sigma, variability %.5lg"%tuple(source_stats)))

                source_results.append([name]+source_stats)
                

            setattr(self,fn.replace(".fits",""),da.DataFile(fn))
            
            if rmf_fn is not None:
                setattr(self,rmf_fn.replace(".fits",""),da.DataFile(rmf_fn))

            if arf_fn is not None:
                setattr(self,arf_fn.replace(".fits",""),da.DataFile(arf_fn))

            self.extracted_sources.append([name,fn.replace(".fits",""),rmf_fn.replace(".fits",""),arf_fn.replace(".fits","")])

        self.source_results=source_results

        srf=open("source_summary.txt","w")
        for sr in source_results:
            srf.write(sr[0].replace(" ","_")+" "+" ".join(["%.5lg"%s for s in sr[1:]])+"\n")

        self.source_summary=da.DataFile("source_summary.txt")
            
 #       for l in allsource_summary:
#            print(l[0] #,l[1].shape)


import dataanalysis.callback

dataanalysis.callback.default_callback_filter.set_callback_accepted_classes([ISGRISpectraSum])
