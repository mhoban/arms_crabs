SHELL := bash
DOCX_FILES := output/$(strip $(patsubst %.Rmd, %.docx, $(wildcard *.Rmd)))
TABLES = $(wildcard data/tables/*.csv)
REFDOC = resources/ref.docx
CITES = references/citations.bib
CITESTYLE = references/style.csl

define OSASCRIPT
tell application "Microsoft Word"
	try
    activate
		close document "$(notdir $(DOCX_FILES))" saving no
    activate
  on error
	end try
end tell
endef

all: docx
force: clean docx

export OSASCRIPT
force-open: clean docx
	@echo "$$OSASCRIPT" | osascript
	@open $(DOCX_FILES)
	
open: docx
	@echo "$$OSASCRIPT" | osascript
	@open $(DOCX_FILES)

docx: $(DOCX_FILES) 

output/%.docx: %.Rmd $(TABLES) $(REFDOC) $(CITES) $(CITESTYLE)
	@echo building $@
	@R --slave -e 'rmarkdown::render("$<",output_file="$@")'

.PHONY: clean
clean:
	@echo cleaning up...
	@$(RM) -f cache/*
	@$(RM) -f $(DOCX_FILES)
	@$(RM) -f output/figures/*.svg