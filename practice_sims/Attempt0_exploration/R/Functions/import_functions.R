#function imports .arp files and converts them to .gen files 
import_arp2gen_files = function(mypath, mypattern) {
  temp_list_1 = list.files(mypath, mypattern, ...)
  temp_list_2 = list(length = length(temp_list_1))
  for(i in 1:length(temp_list_1)){temp_list_2[[i]]=arp2gen(temp_list_1[i])}
  temp_list_2
}
#to call function...
#import_arp2gen_files("pathway\\folders", ".arp$")


#function that converts .gen files to genind objects
import_gen2genind_objects = function(mypath, mypattern) {
  temp_list_3 = list.files(mypath, mypattern)
  temp_list_4 = list(length = length(temp_list_3))
  for(j in 1:length(temp_list_3)){temp_list_4[[j]]=read.genepop(temp_list_3[j], ncode=3)}
  temp_list_4
}
##to call function...
#import_gen2genind_objects("pathway\\folders", ".gen$)