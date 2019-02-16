CPP=g++
LDFLAGS= -I/inc/ -lGL -lGLEW -lglfw
OBJDIR=build/
# Put the dependencies files here
depsname= Window Shader Mesh Camera Texture
foldersrc=src/

# Preparing full name of files & objs
depsobjpref=$(addprefix $(OBJDIR), $(depsname))
depsfilesname=$(addsuffix .cpp, $(depsname))

depsobj=$(addsuffix .o, $(depsobjpref))
depsfiles=$(addprefix $(foldersrc), $(depsfilesname))


TARGETNAME=clouds
TARGET=$(addprefix $(OBJDIR), $(TARGETNAME))
#Window.o: Window.cpp
#	$(CPP) -c $(LDFLAGS) $< -o $@
$(TARGETNAME): $(depsobj)
	$(CPP) src/main.cpp -o $(TARGET) $^ $(LDFLAGS)
.PHONY: $(TARGETNAME)

build/%.o: src/%.cpp
	$(CPP) -c $(LDFLAGS) $< -o $@
#$(depsobj): $(depsfiles)
#	$(CPP) -c $(LDFLAGS) $< -o $@

clean:
	rm -f $(depsobj) $(TARGET)
.PHONY: clean

run:
	./$(TARGET)
.PHONY: run

# Just for testcases
print:
	echo $(TARGET)
.PHONY: print