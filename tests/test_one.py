import sys
import os

__this_dir__ = os.path.join(os.path.abspath(os.path.dirname(__file__)))

import glob

import astropy.io.fits as fits

import ddosa
import findic
import ddosa11
import process_isgri_spectra

# class SpectraBins(ddosa.SpectraBins):
#     def main(self):
#         self.binrmf="/home/isdc/savchenk/osa11_deployment/deployment/ic/ic/ibis/rsp/rmf_62bands.fits"
#         e=fits.open(self.binrmf)[3].data
#         self.bins=list(zip(e['E_MIN'],e['E_MAX']))
#         self.binrmfext=self.binrmf+'[3]'

# class ISGRIResponse(ddosa.ISGRIResponse):
#     path="/home/isdc/savchenk/osa11_deployment/deployment/ic/ic/ibis/rsp/rmf_62bands.fits"


def test_sum():
    #scw = "066500220010.001"
    scw = "198700220010.001"

    fa=process_isgri_spectra.ISGRISpectraSum(assume=[
                            process_isgri_spectra.ScWSpectraList(
                                input_scwlist=ddosa.IDScWList(
                                    use_scwid_list=[scw],
                                    use_version="v1"+scw))])
    fa.cached = False    
    fa.get()

    # print(fa.spectra)

    print(fa.extracted_sources)


def test_one():
    fa=process_isgri_spectra.ProcessSpectra(assume=[
                            ddosa.ScWData(input_scwid="066500220010.001"),
                        ])
    fa.get()
