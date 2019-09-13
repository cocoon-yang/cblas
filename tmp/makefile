# Target   
PROGRAM = cblas.dll  
  
inc = "../inc/"  
src = "../src/"    
INCLUDEDIRS =  /I $(inc)    
 
# Flags  
CPPOPT = $(INCLUDEDIRS) /w /EHsc /D_CRT_SECURE_NO_DEPRECATE  
  
# Compiler   
cc = cl   
CFLAGS =   
   
# The CPPSOURCES macro contains a list of source files.
CPPSOURCES =  cblas.cc      
  
# The CPPOBJECTS macro converts the CPPSOURCES macro into a list
# of object files.         
CPPOBJECTS = $(CPPSOURCES:.cc=.obj)   
    
all: $(PROGRAM)  
  
$(PROGRAM): $(CPPOBJECTS)  
    link.exe /out:$(PROGRAM)  /DLL  $(CPPOBJECTS)    
      
cblas.obj:   
    $(cc) $(CPPOPT) /c ../src/cblas.cc       
               
clean:    
    del $(CPPOBJECTS) $(PROGRAM)  
