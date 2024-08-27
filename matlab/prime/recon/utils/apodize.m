function [ res ] = apodize( img, apodization_para )
%APODIZE Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
    apodization_para=0.15;
end

    kcomb = mrir_fDFT_freqencode(mrir_fDFT_phasencode(img));

    kapodize = mrir_filter_raw_apodize_1d( mrir_filter_raw_apodize_1d(kcomb, 1, apodization_para),  2, apodization_para) ;

    res = mrir_iDFT_freqencode(mrir_iDFT_phasencode(kapodize));
    
end

