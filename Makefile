SHELL := bash
DOCX_FILES := output/$(strip $(patsubst %.Rmd, %.docx, $(wildcard *.Rmd)))

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

output/%.docx: %.Rmd
	@echo building $@
	@R --slave -e 'rmarkdown::render("$<",output_file="$@")'

.PHONY: clean
clean:
	@echo cleaning up...
	@$(RM) -f $(DOCX_FILES)