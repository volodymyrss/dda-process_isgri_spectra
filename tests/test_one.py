import glob

import astropy.io.fits as fits

import ddosa
import process_isgri_spectra

class SpectraBins(ddosa.SpectraBins):
    def main(self):
        self.binrmf="/home/isdc/savchenk/osa11_deployment/deployment/ic/ic/ibis/rsp/rmf_62bands.fits"
        e=fits.open(self.binrmf)[3].data
        self.bins=zip(e['E_MIN'],e['E_MAX'])
        self.binrmfext=self.binrmf+'[3]'

class ISGRIResponse(ddosa.ISGRIResponse):
    path="/home/isdc/savchenk/osa11_deployment/deployment/ic/ic/ibis/rsp/rmf_62bands.fits"


def test_sum():
    fa=process_isgri_spectra.ISGRISpectraSum(assume=[
                            process_isgri_spectra.ScWSpectraList(input_scwlist=ddosa.IDScWList(use_scwid_list=["066500220010.001"],use_version="v1066500220010")),
                        ])
    fa.get()

    print fa.spectra

    print fa.extracted_sources


def test_one():
    fa=process_isgri_spectra.ProcessSpectra(assume=[
                            ddosa.ScWData(input_scwid="066500220010.001"),
                        ])
    fa.get()
