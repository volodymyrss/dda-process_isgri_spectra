

import ddosa 
from astropy.io import fits 

try:
    from pscolors import render
except:
    from bcolors import render

import datetime
import subprocess
import os
import shutil
import pilton


try:
    from dataanalysis import core as da
except ImportError:
    import dataanalysis as da

from numpy import *
from collections import defaultdict

import useresponse

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


class ISGRISpectrumPack(ddosa.DataAnalysis):
    input_spectra=ddosa.ii_spectra_extract
    input_response=useresponse.RebinResponse
          
    cached=True

    def main(self):
        if hasattr(self.input_spectra,'empty_results'):
            print("skipping")
            self.empty_results = True
            return

    
        self.spectra_spectrum = ddosa.localized_DataFile(self.input_spectra.spectrum)
        self.rmf = ddosa.localized_DataFile(self.input_response.rmf)

class ProcessSpectra(ddosa.DataAnalysis):
    input_spectra=ddosa.ii_spectra_extract
    input_response=useresponse.RebinResponse
    #input_response=useresponse.FindResponse
    #input_response=ddosa.ISGRIResponse


    def main(self):
        f=fits.open(self.input_spectra.spectrum.path)


        for hdu in f[2:]:
            sname=hdu.header['NAME']
            fn="isgri_spectrum_%s.fits"%sname.replace(" ","_")

            #hdu.header['ANCRFILE']=self.input_arf.arf_path
            hdu.header['RESPFILE']=self.input_response.rmf_path()
            if hdu.header['RESPFILE'][-1]=="&":
                hdu.header['RESPFILE']=hdu.header['RESPFILE'][:-1]

            if not os.path.exists(hdu.header['RESPFILE']):
                raise Exception("no rsp file")

            print("source:",sname,"to",fn)
            print("response",hdu.header['RESPFILE'])
            print(hdu.header['ANCRFILE'])
            hdu.writeto(fn,clobber=True)
            setattr(self,fn,da.DataFile(fn))



class ScWSpectraList(ddosa.DataAnalysis):
    input_scwlist=ddosa.RevScWList
    copy_cached_input=False
    input_specsummary=ddosa.SpectraProcessingSummary
    input_rebinresponsesummary=useresponse.RebinResponseProcessingSummary

    allow_alias=True

    version="allthem"
    
    maxspec=None

    def main(self):
        self.spectra=[ISGRISpectrumPack(assume=scw) for scw in self.input_scwlist.scwlistdata]
        #self.spectra=[[
                    #    ddosa.ii_spectra_extract(assume=scw),
                    #    useresponse.RebinResponse(assume=scw),
                    #    ddosa.ISGRIResponse(assume=scw)
                    #    ] for scw in self.input_scwlist.scwlistdata]

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

        e1,e2=list(map(array,list(zip(*self.input_bins.bins))))
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

    spectra=None

    cached=True

    version="v5.8.2"

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

#        for spectrum, rmf, arf in choice:
#            print("processing",spectrum,rmf,arf)
        for pack in choice:
            if hasattr(pack,'empty_results'):
                print("skipping", pack)
                continue
            
            fn=pack.spectra_spectrum.get_path()
            print("%i/%i"%(i_spec,len(choice)))
            tc=time.time()
            print("seconds per spectrum:",(tc-t0)/i_spec,"will be ready in %.5lg seconds"%((len(choice)-i_spec)*(tc-t0)/i_spec))
            i_spec+=1
            print("spectrum from",fn)

            f=fits.open(fn)


            t1,t2=f[1].header['TSTART'],f[1].header['TSTOP']
            print(t1,t2)

            try:
                print("found extensions",len(list(f)))
            except Exception in e:
                print("failed to see extensions",e)
                continue

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
                        spectra[name] = [rate, err**2, exposure, e.copy(), defaultdict(int), defaultdict(int),ontime,telapse,tstart,tstop,[revol],exp_src]
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

                    if hasattr(pack.arf, 'arf_path'):
                        arf_path = pack.arf.get_path()
                    else:
                        arf_path = None

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

        eb1,eb2=list(map(array,list(zip(*self.input_response.bins))))

        self.spectra=spectra

        source_results=[]
        self.extracted_sources=[]

        for name,spectrum in list(spectra.items()):
            source_short_name=name.strip().replace(" ","_").replace("/","_")

            assert(len(list(spectrum[4].keys()))==1)

            if list(spectrum[4].keys())[0] is not None:
                arf_fn="arf_sum_%s.fits"%source_short_name
                fits.open(list(spectrum[4].keys())[0]).writeto(arf_fn,clobber=True)
            else:
                arf_fn="arf_sum_%s.fits"%source_short_name
                dc=pilton.heatool("dal_create")
                dc['obj_name']=arf_fn
                dc['template']="ISGR-ARF.-RSP.tpl"

                ddosa.remove_withtemplate(arf_fn+"(ISGR-ARF.-RSP.tpl)")
                dc.run()

                _rmf=fits.open(list(spectrum[5].keys())[0])

                _arf=fits.open(arf_fn)
                _arf[1].data=zeros(len(_rmf['ISGR-RMF.-RSP'].data['ENERG_LO']),dtype=_arf[1].data.dtype)

                _arf[1].data['ENERG_LO']=_rmf['ISGR-RMF.-RSP'].data['ENERG_LO']
                _arf[1].data['ENERG_HI']=_rmf['ISGR-RMF.-RSP'].data['ENERG_HI']
                _arf[1].data['SPECRESP']=1.
                _arf.writeto(arf_fn,overwrite=True)

                #fits.open(arf_fn)
                #arf_fn=None
            
            print("response keys",len(list(spectrum[5].keys())),list(spectrum[5].keys()))
            spectrum[5]=dict(list(spectrum[5].items())[:1])

            assert(len(list(spectrum[5].keys()))==1)
            rmf_fn="rmf_sum_%s.fits"%source_short_name

            print("writing:",list(spectrum[5].keys())[0],"to",rmf_fn)
            fits.open(list(spectrum[5].keys())[0]).writeto(rmf_fn,clobber=True)
            

            try:
                spectrum[3].data['RATE'][:] = spectrum[0][:]
                spectrum[3].data['STAT_ERR'][:] = (spectrum[1]**0.5)[:]
                #spectrum[3].data['RATE'][:],spectrum[3].data['STAT_ERR'][:] = self.input_efficiency.correct(spectrum[0][:],(spectrum[1]**0.5)[:])
            except Exception as e:
                print("problem with", spectrum[0], spectrum[1], spectrum[3], e)
                raise

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


            _rmf=fits.open(list(spectrum[5].keys())[0])
            spectrum[3].data['QUALITY'][_rmf[1].data['E_MIN']<30]=3

            fn="isgri_sum_%s.fits"%source_short_name
            spectrum[3].writeto(fn, clobber=True, checksum=True)
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

                lc_t1,lc_t2,lc_f,lc_fe=list(map(array,list(zip(*lc))))

                if self.save_lc:
                    savetxt("%s_%.5lg_%.5lg.txt"%(source_short_name,erange[0],erange[1]),column_stack((lc_t1,lc_t2,lc_f,lc_fe)))

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

            self.extracted_sources.append([name,fn.replace(".fits",""),rmf_fn.replace(".fits",""),arf_fn.replace(".fits","") if arf_fn is not None else ""])

        self.source_results=source_results

        srf=open("source_summary.txt","w")
        for sr in source_results:
            srf.write(sr[0].replace(" ","_")+" "+" ".join(["%.5lg"%s for s in sr[1:]])+"\n")


        #heaspa.PHA(total_spectrum_summed,total_spectrum_summed**0.5,exposure=total_exposure).write(total_fn)
            
 #       for l in allsource_summary:
#            print(l[0] #,l[1].shape)



import dataanalysis.callback

dataanalysis.callback.default_callback_filter.set_callback_accepted_classes([ddosa.mosaic_ii_skyimage, ddosa.ii_skyimage, ddosa.BinEventsImage, ddosa.ibis_gti, ddosa.ibis_dead, ddosa.ISGRIEvents, ddosa.ii_spectra_extract, ddosa.BinEventsSpectra, ddosa.ii_lc_extract, ddosa.BinEventsLC, ISGRISpectraSum])


