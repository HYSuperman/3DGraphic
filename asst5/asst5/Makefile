BASE = asst5

all: $(BASE)

OS := $(shell uname -s)

ifeq ($(OS), Linux) # Science Center Linux Boxes
  CPPFLAGS = -I/home/l/i/lib175/usr/glew/include
  LDFLAGS += -L/home/l/i/lib175/usr/glew/lib -L/usr/X11R6/lib
  LIBS += -lGL -lGLU -lglut
endif

ifeq ($(OS), Darwin) # Assume OS X
  CPPFLAGS += -D__MAC__
  CPPFLAGS += -I/usr/include -I/usr/local/include
  LDFLAGS += -framework GLUT -framework OpenGL
  # suppress the warnings about deprecated glu/glut functions on Mavericks
  CXXFLAGS += -Wno-deprecated
endif

ifdef OPT 
  #turn on optimization
  CXXFLAGS += -O2
else 
  #turn on debugging
  CXXFLAGS += -g
endif

CXX = g++ 

OBJ = $(BASE).o ppm.o glsupport.o scenegraph.o picker.o

$(BASE): $(OBJ)
#	$(LINK.cpp) -o $@ $^ $(LIBS) -lGLEW
#	Modified using my local library path
	$(LINK.cpp) -o $@ $^ $(LIBS) /usr/local/lib/libGLEW.a

clean:
	rm -f $(OBJ) $(BASE)
