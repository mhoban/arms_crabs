local utils = require 'pandoc.utils'
local run_json_filter = utils.run_json_filter

-- This works on newer Pandoc versions but doesn't on pandoc 2.2.3.2
local function run_citeproc (doc)
  if PANDOC_VERSION >= '2.19.1' then
    return pandoc.utils.citeproc(doc)
  elseif PANDOC_VERSION >= '2.11' then
    local args = {'--from=json', '--to=json', '--citeproc'}
    return run_json_filter(doc, 'pandoc', args)
  else
    return run_json_filter(doc, 'pandoc-citeproc', {FORMAT, '-q'})
  end
end

-- local function run_citeproc (doc)
  -- return run_json_filter(doc, 'pandoc-citeproc')
-- end

--- Filter to the references div and bibliography header added by
--- pandoc-citeproc.
local remove_pandoc_citeproc_results = {
  Header = function (header)
    return header.identifier == 'bibliography'
      and {}
      or nil
  end,
  Div = function (div)
    return div.identifier == 'refs'
      and {}
      or nil
  end
}

-- stackoverflow
function table.contains(table, element)
  for _, value in pairs(table) do
    if value == element then
      return true
    end
  end
  return false
end

function create_bibliographies (doc)
  local blocks = {}
  local new_blocks = {}
  for block_id,block_data in pairs(doc.blocks) do
    if block_data.attr and block_data.attr.classes and table.contains(block_data.attr.classes, "multi-refs") then
      local tmp_doc = pandoc.Pandoc(new_blocks, doc.meta)
      local new_doc = run_citeproc(tmp_doc)
      for _, block_to_add in pairs(new_doc.blocks) do
        blocks[#blocks+1] = block_to_add
      end
      new_blocks = {}
    else
      new_blocks[#new_blocks+1] = block_data
    end
  end
  for _, new_block in pairs(new_blocks) do
    blocks[#blocks+1] = new_block
  end
  return pandoc.Pandoc(blocks, doc.meta)
end

return {
  -- remove result of previous pandoc-citeproc run (for backwards
  -- compatibility)
  remove_pandoc_citeproc_results,
  {Pandoc = create_bibliographies},
}
