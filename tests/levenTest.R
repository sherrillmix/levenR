##File name: levenTest.R
##Creation date: Oct 01, 2012
##Last modified: Mon Oct 01, 2012  10:00AM
##Created by: scott
##Summary: Test leven functions 
library(leven)
if(leven('AAA','AAA')!=0)stop(simpleError('Equal strings not equal'))
if(leven('AT','A')!=1)stop(simpleError('Difference between AT and A not 1'))
if(leven('A','AT')!=1)stop(simpleError('Difference between A and AT not 1'))
if(leven('A','T')!=1)stop(simpleError('Difference between A and T not 1'))
if(leven('asdasdasedaeswaseaseasdasdfdfsdfgdfgsdfsdfaeruiurtyuuiyuyiuvdyiuyiuyiuyiuyiuyiuyiuyasduifysaiyufasiyudfaisyufiasyufisayufiasyufisayfisayufsaiyufisaydfiasyufiasydfiasyfisaydfiasyisayfisyisaiuyyiauyriutryiuayiuyiuyiuyiuyiuyiuyiuyiuadfuasduyasfiuysadufyaisydfiuaysfiyiuyiuaydfiasydfiasyayyig','asdasdasedaeswaseaseasdasdfdfsdfgdfgsdfsdfaeruiurtyuuiyuyiuvdyiuyiuyiuyiuyiuyiuyiuyasduifysaiyufasiyudfaisyufiasyufisayufiasyufisayfisayufsaiyufisaydfiasyufiasydfiasyfisaydfiasyisayfisyisaiuyyiauyriutryiuayiuyiuyiuyiuyiuyiuyiuyiuadfuasduyasfiuysadufyaisydfiuaysfiyiuyiuaydfiasydfiasyayyig')!=0)stop(simpleError('Equal longer strings not equal'))
if(all(leven(c('AAAA','AA'),'AAA')!=1))stop(simpleError('Multiple string1 do not produce expected'))
if(all(leven('AAA',c('AAAA','AA'))!=1))stop(simpleError('Multiple string2 do not produce expected'))
if(tryCatch(leven(NULL,'AAA'),error=function(x)"Caught error")!='Caught error')stop(simpleError('Did not catch NULL input'))
if(all(unlist(levenAlign('AA','AA'))!=c('AA','AA')))stop(simpleError('Aligning AA and AA does not produce AA and AA'))
if(all(unlist(levenAlign('AAT','AA'))!=c('AA-','AAT')))stop(simpleError('Aligning AAT and AA does not produce AAT and AA-'))
if(all(unlist(levenAlign('AATC','AAC'))!=c('AA-C','AATC')))stop(simpleError('Aligning AATC and AAC does not produce AATC and AA-C'))
if(all(unlist(levenAlign('AAC','AATC'))!=c('AATC','AA-C')))stop(simpleError('Aligning AAC and AATC does not produce AA-C and AATC'))
if(all(unlist(levenAlign('TAAAAAATC','AAATC'))!=c('-AAA---TC','TAAAAAATC')))stop(simpleError('Aligning TAAAAAATC and AAATC does not produce -AAA---TC and TAAAAAATC'))
if(all(unlist(levenAlign(c('TAA','AAT'),'AA'))!=c('-AA-','TAA-','-AAT')))stop(simpleError('Aligning TAA,AAT and AA does not produce TAA-,-AAT and -AA-'))
