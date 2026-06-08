USE_MKL ?= 0

ifeq ($(USE_MKL), 1)
    CXX      := icpx
    CXXFLAGS := -Wall -c -qopenmp -O2 -DUSE_MKL
    LDFLAGS  := -Wall -qmkl=parallel
else
    CXX      := g++
    CXXFLAGS := -Wall -c -fopenmp -O2
    LDFLAGS  := -Wall -fopenmp
endif

CXXVERSION := -std=c++17
OBJDIR := obj
BINDIR := ./
DISTDIR := dist
EXEFILE := mc_sim

SRCDIR := src
verz := $(foreach dir,$(SRCDIR),$(wildcard $(dir)/*))
cppfiles := $(foreach dir,$(verz),$(wildcard $(dir)/*.cpp))

objects := $(subst $(SRCDIR),$(OBJDIR),$(cppfiles:.cpp=.o))
deps := $(objects:.o=.d)

.PHONY: all clean dist

all: $(BINDIR)/$(EXEFILE)

-include $(deps)

$(BINDIR)/$(EXEFILE): $(objects)
	mkdir -p $(BINDIR)
	$(CXX) $(CXXVERSION) $(LDFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo compiling  $< ...
	$(CXX) $(CXXVERSION) $(CXXFLAGS) $< -o $@

$(OBJDIR)/%.d: $(SRCDIR)/%.cpp
	mkdir -p $(@D)
	$(CXX) $(CXXVERSION) -MM -MT "$@ $(patsubst %.d,%.o,$@)" -MF $@ $<

$(DISTDIR)/gepackt.tar.gz: $(BINDIR)/$(EXEFILE)
	mkdir -p $(DISTDIR)
	tar cvf $(DISTDIR)/gepackt.tar $(BINDIR)/$(EXEFILE)
	gzip $(DISTDIR)/gepackt.tar

dist: $(DISTDIR)/gepackt.tar.gz

clean:
	$(RM) -rf $(OBJDIR)
	$(RM) $(EXEFILE)

ultraclean: clean
	$(RM) -rf $(DISTDIR)
