# Target     
PROGRAM = lu.exe    
    
LIBDIRS = "../lib/" 
inc = "../include/"    
src = "../src/"      
INCLUDEDIRS =  /I $(inc)  /I $(LIBDIRS)    
   
# Flags    
CPPOPT = $(INCLUDEDIRS) /w /EHsc /D_CRT_SECURE_NO_DEPRECATE /D ADD_ /D HAVE_LAPACK_CONFIG_H \
	/D LAPACK_COMPLEX_STRUCTURE /D WIN32 /D NDEBUG /D _CONSOLE    
    
 
LIBS =  ../lib/liblapacke.lib  ../lib/liblapack.lib  ../lib/libblas.lib ../lib/cblas.lib

# Compiler     
cc = cl     
CFLAGS =     
     
LIBFLAGS =  /LIBPATH $(LIBDIRS)  

# list of source files     
CPPSOURCES =  lu.cpp  foo.cpp     
    
# expands to list of object files            
CPPOBJECTS = $(CPPSOURCES:.cpp=.obj)     
      
all: $(PROGRAM)    
    
$(PROGRAM): $(CPPOBJECTS)  # $(LIBS)  
    link.exe /out:$(PROGRAM)  $(CPPOBJECTS)   $(LIBS)                     
  
lu.obj:     
     $(cc) $(CPPOPT) /c ../src/lu.cpp          

foo.obj:     
     $(cc) $(CPPOPT) /c ../src/foo.cpp   

clean:      
    del $(CPPOBJECTS) $(PROGRAM)  
