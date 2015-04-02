pro check


;Restore the SAV files
RESTORE, '/home/omiros/SREMDC/process/output/FSVD.sav'
RESTORE, '/home/omiros/SREMDC/process/output/FSVD_old.sav'

check1=fluxes_info.FPDO[1,10]-fluxes_info_old.FPDO[1,10]
check2=fluxes_info.FPDO[5,10]-fluxes_info_old.FPDO[5,10]
check3=fluxes_info.FEDO[1,10]-fluxes_info_old.FEDO[1,10]

print, check1, check2, fluxes_info.FPDO[1,10], fluxes_info.FEDO[1,10]


stop
end


