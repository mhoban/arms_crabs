SHELL := bash
DOCX_FILES := output/$(strip $(patsubst %.Rmd, %.docx, $(wildcard *.Rmd)))
MD_FILES := output/$(strip $(patsubst %.Rmd, %.Rmd.md, $(wildcard *.Rmd)))
R_FILES := crabs.R lib/CAPdiscrim.R collapse_taxonomy.R
TABLES = $(wildcard data/tables/*.csv)
REFDOC = resources/ref.docx
CITES = references/citations.json
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

output/%.docx: %.Rmd $(TABLES) $(REFDOC) $(CITES) $(CITESTYLE) $(R_FILES)
	@echo building $@
	@R --slave -e 'knitr::knit("$<","output/$<.md")'
	@pandoc +RTS -K512m -RTS --filter=pandoc-crossref \
		--citeproc  output/$<.md  --to docx \
		--from markdown+autolink_bare_uris+tex_math_single_backslash \
		--output $@ --lua-filter resources/pagebreak.lua \
		--highlight-style tango --toc-depth=2 \
		--lua-filter=resources/scholarly-metadata.lua \
		--lua-filter=resources/author-info-blocks.lua \
		--lua-filter=resources/abstract-section.lua \
		--reference-doc=resources/ref.docx \
		--lua-filter=resources/multirefs.lua

.PHONY: clean
clean:
	@echo cleaning up...
	@$(RM) -f cache/*
	@$(RM) -f $(DOCX_FILES)
	@$(RM) -f $(MD_FILES)
	@$(RM) -f output/figures/*.svg