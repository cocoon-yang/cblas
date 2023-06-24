# Target   
PROGRAM = cblas.dll  
  
inc = "../inc/"  
src = "../src/"    
INCLUDEDIRS =  /I $(inc)    
 
# Flags  
CPPOPT = $(INCLUDEDIRS) /w /EHsc /O1 /D_CRT_SECURE_NO_DEPRECATE  
  
# Compiler   
cc = cl   
CFLAGS =   
   
# The CPPSOURCES macro contains a list of source files.
CPPSOURCES =  cblas1_d.cc    cblas2_d.cc    cblas3_d.cc   lapack_d.cc    assistant.cc
  
# The CPPOBJECTS macro converts the CPPSOURCES macro into a list
# of object files.         
CPPOBJECTS = $(CPPSOURCES:.cc=.obj)   
    
all: $(PROGRAM)  
  
$(PROGRAM): $(CPPOBJECTS)  
    link.exe /out:$(PROGRAM)  /DLL  $(CPPOBJECTS)    
      
assistant.obj:   
    $(cc) $(CPPOPT) /c ../src/assistant.cc     

cblas1_d.obj:   
    $(cc) $(CPPOPT) /c ../src/cblas1_d.cc       
          
cblas2_d.obj:   
    $(cc) $(CPPOPT) /c ../src/cblas2_d.cc     

cblas3_d.obj:   
    $(cc) $(CPPOPT) /c ../src/cblas3_d.cc     
       
lapack_d.obj:   
    $(cc) $(CPPOPT) /c ../src/lapack_d.cc     
    
clean:    
    del $(CPPOBJECTS) $(PROGRAM)  

update:
    copy  $(PROGRAM)  ..\test\$(PROGRAM)  