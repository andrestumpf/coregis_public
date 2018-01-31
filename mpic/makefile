# Yet another makefile ...

TARGET   = MPIC
CC       = gcc 
#CC       = icc 

GCC_FLAGS = -std=gnu99 -O3 -fopenmp -ffast-math -msse2 -Wall -I. -I/usr/include/x86_64-linux-gnu -I/usr/include/gdal
ICC_FLAGS = -std=gnu99 -O3 -fopenmp -no-prec-div -ansi-alias -xCORE-AVX2 -Wall -I. -I/usr/include/x86_64-linux-gnu -I/usr/include/gdal

LINKER   = gcc -o
#LINKER   = icc -o
LFLAGS   = -Wall -fopenmp -lm -lpcre -ltiff -lgdal
CFLAGS   = $(GCC_FLAGS)

SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

SOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

all : $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(OBJECTS)
	@$(LINKER) $@ $(OBJECTS) $(LFLAGS) 
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c 
	@$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

clean :
	@rm -f $(OBJDIR)/*.o

proper : clean
	@rm -f $(BINDIR)/$(TARGET)
