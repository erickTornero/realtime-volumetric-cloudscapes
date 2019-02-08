CPP=g++
LDFLAGS= -I/inc/ -lGL -lGLEW -lglfw
OBJDIR=build/
depsname= Window
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

$(depsobj): $(depsfiles)
	$(CPP) -c $(LDFLAGS) $< -o $@


.PHONY: clouds

clean:
	rm -f $(depsobj) $(TARGET)
.PHONY: clean

run:
	./$(TARGET)

print:
	echo $(TARGET)