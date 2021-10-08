import ddosa 
from astropy.io import fits 
import time
from pathlib import Path

try:
    from pscolors import render
except:
    from bcolors import render

import datetime
import subprocess
import os

import pilton

from astropy.time import Time

try:
    from dataanalysis import core as da
except ImportError:
    import dataanalysis as da

from numpy import *
import numpy as np
from collections import defaultdict

import useresponse


try:
    import crab
except:
    pass

try:
    import matplotlib as mpl
    mpl.use('agg')
    import matplotlib.pylab as plt
except:
    print("no matplotlib: some functionality will not be available")


try:
    import ogip
    import ogip.spec
except ImportError:
    print("no ogip: maybe it's ok")

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
        self.spectra=[[
                       ddosa.ii_spectra_extract(assume=scw),
                       useresponse.RebinResponse(assume=scw),
                       ddosa.ISGRIResponse(assume=scw)
                       ] for scw in self.input_scwlist.scwlistdata]

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


def merge_rmfs(rmfs: dict):
    merged_rmf_f = None
    total_exposure = 0

    for rmf_fn, rmf_exposure in rmfs.items():
        if merged_rmf_f is None:
            merged_rmf_f = fits.open(rmf_fn)
            merged_rmf_f['ISGR-RMF.-RSP'].data['MATRIX'] *= rmf_exposure
            total_exposure += rmf_exposure
        else:
            rmf_f = fits.open(rmf_fn)
            merged_rmf_f['ISGR-RMF.-RSP'].data['MATRIX'] + \
                rmf_f['ISGR-RMF.-RSP'].data['MATRIX'] * rmf_exposure
            total_exposure += rmf_exposure

    merged_rmf_f['ISGR-RMF.-RSP'].data['MATRIX'] /= total_exposure

    print("merged rmf total exposure", total_exposure)

    return merged_rmf_f, total_exposure

class ISGRISpectraSum(ddosa.DataAnalysis):
    version="v5.8.6"

    input_spectralist=ScWSpectraList
    input_response=ddosa.SpectraBins
    input_efficiency=SpectrumEfficiencyCorrection

    input_icroot=ddosa.ICRoot

    copy_cached_input=False

    verify_background_lines=True
    collect_ic=True
    adapt_zero_offset=False

    spectra=None

    cached=True

    sources=['Crab']

    extract_all=True
    save_lc=True

    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.extract_all:
            v+=".extractall"
        else:
            v+".extract_"+("_".join(self.sources))

        if self.adapt_zero_offset:
            v += ".offsetmorph"

        return v

    def main(self):
        spectra={}

        choice=self.input_spectralist.spectra

        allsource_summary=[]


        sig = lambda x, y: (((x/y)[~np.isnan(y) & (y!=0) & ~np.isinf(y) & ~np.isinf(x) & ~np.isnan(x)])**2).sum()**0.5

        t0=time.time()
        i_spec=1

        used_spectra = []

        for spectrum, rmf, arf in choice:
            print("\n\n\033[31mprocessing: \n" + \
                  "          \033[32m" + str(spectrum) + "\n"\
                  "          \033[32m" + str(rmf) + "\n"\
                  "          \033[32m" + str(arf) + "\033[0m")

            if hasattr(spectrum,'empty_results'):
                print("skipping",spectrum)
                continue

            if not hasattr(spectrum,'spectrum'):
                print("skipping",spectrum)
                continue
            
            fn = spectrum.spectrum.get_path()
            rmf_fn = rmf.rmf.get_path()

            print("%i/%i"%(i_spec,len(choice)))
            tc=time.time()
            print("seconds per spectrum:",(tc-t0)/i_spec,"will be ready in %.5lg seconds"%((len(choice)-i_spec)*(tc-t0)/i_spec))
            i_spec+=1

            if fn in used_spectra:
                raise Exception("repeated spectrum! BIG problem with environment!")

            used_spectra.append(fn)
            print("\033[33mspectrum from",fn,"\033[0m")
            print("\033[33mrmf from",rmf_fn,"\033[0m")

            f=fits.open(fn)
            f_rmf=fits.open(rmf_fn)

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
                        # new
                        spectra[name] = [rate, err**2, exposure, e.copy(), defaultdict(int), defaultdict(int),
                                         ontime, telapse, tstart, tstop, [revol], exp_src]
                    else:
                        # add
                        err[np.isnan(err) | (err==0)]=np.inf
                        spectra[name][1][np.isnan(spectra[name][1]) | (spectra[name][1]==0)]=np.inf
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

                    if hasattr(arf,'arf_path'):
                        arf_path=arf.arf_path
                    else:
                        arf_path=None

                    rmf_path=rmf.rmf.get_path()

                    spectra[name][4][arf_path]+=exposure
                    spectra[name][5][rmf_path]+=exposure

                    print(f"\033[31mcurrent exposure for {rmf_path} {spectra[name][5][rmf_path]}\033[0m")
            
                    print(render("{BLUE}%.20s{/}"%name),"%.4lg sigma in %.5lg ks"%(sig(rate,err),exposure/1e3),"total %.4lg in %.5lg ks"%(sig(spectra[name][0],spectra[name][1]), spectra[name][2]/1e3))

            f.close()

            try:
                print(get_open_fds())
            except Exception as e:
                print("unable to check open fds")

        eb1, eb2 = list(map(array,list(zip(*self.input_response.bins))))

        # now write        
        self.spectra=spectra

        source_results=[]
        self.extracted_sources=[]

        self.multi_spec_plot_init()

        for name, spectrum in list(spectra.items()):
            source_short_name=name.strip().replace(" ","_").replace("/","_")

            assert(len(list(spectrum[4].keys()))==1)

            arf_fn = f"arf_sum_{source_short_name}.fits"
            if list(spectrum[4].keys())[0] is not None:
                fits.open(list(spectrum[4].keys())[0]).writeto(arf_fn, clobber=True)
            else:
                self.write_unitary_arf(rmf_fn, arf_fn)
            
            print("response keys",len(list(spectrum[5].keys())),list(spectrum[5].keys()))

            merged_rmf_f, merged_rmf_exposure = merge_rmfs(spectrum[5])
                        
            rmf_fn="rmf_sum_%s.fits"%source_short_name

            print("writing:", merged_rmf_f, "to", rmf_fn)
            merged_rmf_f.writeto(rmf_fn, clobber=True)

            e_morph, e_morph_meta = self.produce_e_morph(
                fits.open(fn)[1].header['TSTART'],
                fits.open(fn)[1].header['TSTOP']
            )

            self.morph_rmf(rmf_fn, arf_fn, e_morph, e_morph_meta)

            

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

            try:
                ebds_ext = _rmf['ISGR-EBDS-MOD']
            except:
                ebds_ext = _rmf['EBOUNDS']

            t = Time(spectrum[3].header['TSTOP']+51544, format='mjd')
            emin = 15 + (t.byear - 2002)/(2020 - 2002) * (30 - 15)

            print("\033[31mselected emin: ", emin, "\033[0m")

            spectrum[3].data['QUALITY'][ebds_ext.data['E_MIN'] < emin]=3

            fn="isgri_sum_%s.fits"%source_short_name
            spectrum[3].writeto(fn, clobber=True, checksum=True)
            print("writing",fn)


            select_range=lambda x,a,b:((eb1>a) & (eb2<b) & ~np.isnan(x) & ~np.isinf(x))

            print(render("{RED}total{/} for {BLUE}%s{/}"%name))

            all_spectra=[l for l in allsource_summary if l[0]==name]
            
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

            self.extracted_sources.append([name,
                                           fn.replace(".fits",""),
                                           rmf_fn.replace(".fits",""),
                                           arf_fn.replace(".fits","") if arf_fn is not None else ""])

            
            if name == 'Background':
                self.asses_background_lines(fn, rmf_fn)

            self.multi_spec_plot_add_spectrum(fn, rmf_fn, e_morph)

        try:
            self.multi_spec_plot_finalize()
        except Exception as e:
            print("unable to plot:", e)

        self.source_results=source_results

        srf=open("source_summary.txt","w")
        for sr in source_results:
            srf.write(sr[0].replace(" ","_")+" "+" ".join(["%.5lg"%s for s in sr[1:]])+"\n")


    def collect_ic_info(self, coverage):
        pass

    def asses_background_lines(self, fn, rmf_fn):
        #TODO: overplot offset and zero-offset morph here

        spec = ogip.spec.PHAI.from_file_name(fn)
                
        # what's  wrong with it?.. compressed?
        #rmf = ogip.spec.RMF.from_file_name(rmf_fn)

        ebds_mod = fits.open(rmf_fn)['ISGR-EBDS-MOD']
        e1 = ebds_mod.data['E_MIN']
        e2 = ebds_mod.data['E_MAX']

        
        t1 = fits.open(fn)[1].header['TSTART']
        t2 = fits.open(fn)[1].header['TSTOP']
        tc = (t2 + t1)/2.

        
        icroot = Path(self.input_icroot.icroot)
                
        lut2idx = fits.open(icroot / 'idx/ic/ISGR-RISE-MOD-IDX.fits')[1].data
        
        p = (icroot / 'ic/ibis' / min(lut2idx, key=lambda r:abs(r['VSTART']-tc))['MEMBER_LOCATION']).resolve()
        lut2 = fits.open(p)[1]

        print(lut2)

        plt.figure()
        plt.imshow(np.stack(lut2.data['CORR']))
        plt.savefig("lut2.png")


        plt.figure(figsize=(10,7))

        plt.plot(
            e1,
            spec._rate/(e2 - e1)
        )
        plt.axvline(59, lw=3, alpha=0.5, c='r')
        plt.axvline(511, lw=3, alpha=0.5, c='r')

        plt.xlim([15, 700])
        #plt.ylim([spec._rate.])
        plt.loglog()
        plt.grid()
        plt.xlabel('keV')
        plt.ylabel('counts/keV')
        plt.savefig("background.png")
        

    def write_unitary_arf(self, rmf_fn, arf_fn):
        dc = pilton.heatool("dal_create")
        dc['obj_name']=arf_fn
        dc['template']="ISGR-ARF.-RSP.tpl"

        ddosa.remove_withtemplate(arf_fn+"(ISGR-ARF.-RSP.tpl)")
        dc.run()

        _rmf=fits.open(rmf_fn)

        _arf=fits.open(arf_fn)
        _arf[1].data=zeros(len(_rmf['ISGR-RMF.-RSP'].data['ENERG_LO']),dtype=_arf[1].data.dtype)

        _arf[1].data['ENERG_LO']=_rmf['ISGR-RMF.-RSP'].data['ENERG_LO']
        _arf[1].data['ENERG_HI']=_rmf['ISGR-RMF.-RSP'].data['ENERG_HI']
        _arf[1].data['SPECRESP']=1.
        _arf.writeto(arf_fn,overwrite=True)


    def produce_e_morph(self, t1, t2):
        if (t2 - t1) > 60:
            raise RuntimeError('too long to patch rmf!')            
            
        from astropy.time import Time

        # need to read from IC, per spectrum, in oda/osa
        # https://gitlab.astro.unige.ch/integral/cc-workflows/cc-global-summary
        def offset_approximation(ijd):
            z1 = 2000
            z2 = 6000
            p1 = -1.5
            p2 = -4
            return (z1 - ijd) * (z2 + (z2-z1) - ijd) / (z1-z2)**2 * (-p2+p1) + p1

        offset_approximation(Time(2015, format='byear').mjd - 51544)

        plt.figure(figsize=(15, 10))

        
        def e_morph_family(E, lebias=4):
            #r = E - E / ( 1 + (E/18)**2 )
            r = E - lebias/(1 + np.exp((E/27)**2))*(1 + np.exp(1))
            return r

        plt.scatter([27, 60, 511], [-4, 0, 0], s=100)

        plt.ylim(-20,20)
        plt.semilogx()

        e = np.logspace(1, 3, 500)

        offset_approx = -offset_approximation((t1+t2)/2.)
        e_morph = lambda en: e_morph_family(en, offset_approx)

        for year in np.linspace(2002, 2022, 10):
            leb = -offset_approximation(Time(year, format='byear').mjd - 51544)
            plt.plot(
                e,
                e_morph_family(e, leb) - e 
            )
            
        plt.plot(
            e,
            e_morph(e) - e,
            lw=5
        )

        plt.grid()
        plt.savefig('ebias-morph.png')
        return e_morph, dict(offset_approx=offset_approx)


    def morph_rmf(self, rmf_fn, arf_fn, e_morph, e_morph_meta):
        import ogip.spec
        import ogip.tools

        rmf = ogip.spec.RMF.from_file_name_osaisgri(rmf_fn) 

        rmf._matrix *= 1.49 # see https://gitlab.astro.unige.ch/integral/cc-workflows/cc-global-summary/ 

        ogip.tools.transform_rmf(
            rmf=rmf,
            arf=ogip.spec.ARF.from_file_name_osa(arf_fn),
            bias_function=e_morph,
            preserved_projection=ogip.tools.crab_ph_cm2_s_kev
        ).to_fits(rmf_fn)
        
        rmf = fits.open(rmf_fn)
        rmf[1].header['MORPH'] = "; ".join({f'{k}: {v}' for k, v in e_morph_meta.items()})
        rmf['EBOUNDS'].header['EXTNAME'] = 'ISGR-EBDS-MOD'
        rmf['MATRIX'].header['EXTNAME'] = 'ISGR-RMF.-RSP'
        rmf.writeto(rmf_fn, overwrite=True)    
                
        
    def multi_spec_plot_init(self):
        self._multi_spec_plot_figure = plt.figure(figsize=(15,10))
        self._max_rate = 0

        
    def multi_spec_plot_add_spectrum(self, fn, rmf_fn, e_morph):
        saved_f = plt.gcf()
        plt.figure(self._multi_spec_plot_figure.number)

        e1 = fits.open(rmf_fn)['ISGR-EBDS-MOD'].data['E_MIN']
        e2 = fits.open(rmf_fn)['ISGR-EBDS-MOD'].data['E_MAX']
        s = fits.open(fn)[1]

        r = s.data['RATE']/(e2-e1)#*e1
        r_morph = s.data['RATE']/(e_morph(e2)-e_morph(e1))#*e1

        x = plt.plot(
            e1,
            r,
        )
        plt.plot(
            e_morph(e1),
            r_morph,
            ls=":",
            c=x[0].get_color(),
            label=f"morphed {s.header['NAME']}"
        )
        self._max_rate = max(np.nanmax(r), self._max_rate)
        
        if s.header['NAME'] == 'Background':
            m_line = (e1>50) & (e1<70)
            
            e1_peak_estim = e1[m_line][np.nanargmax(r[m_line])]
            plt.axvline(e1_peak_estim,
                        label=f'peak estimate: {e1_peak_estim}',
                        lw=5,
                        alpha=0.3
                       )

        plt.figure(saved_f.number)
            

    def multi_spec_plot_finalize(self):
        saved_f = plt.gcf()

        plt.figure(self._multi_spec_plot_figure.number)
        plt.axvline(59, ls=":", lw=3, c='k')
        plt.axvline(511, ls=":", lw=3, c='k')
        plt.axvline(32, ls="--", lw=5, alpha=0.3, c='r')
        plt.axvline(27, ls="--", lw=5, alpha=0.3, c='r')
                
                
        plt.ylim([self._max_rate*0.001, self._max_rate*2])
        plt.loglog()
        plt.legend()

        plt.xlabel('keV')
        plt.ylabel('counts/s/keV')

        plt.savefig("mutli-spec.png")
        plt.figure(saved_f.number)

import dataanalysis.callback

dataanalysis.callback.default_callback_filter.set_callback_accepted_classes([ddosa.mosaic_ii_skyimage, ddosa.ii_skyimage, ddosa.BinEventsImage, ddosa.ibis_gti, ddosa.ibis_dead, ddosa.ISGRIEvents, ddosa.ii_spectra_extract, ddosa.BinEventsSpectra, ddosa.ii_lc_extract, ddosa.BinEventsLC, ISGRISpectraSum])


