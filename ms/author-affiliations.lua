-- Render Quarto's normalized author/affiliation metadata in Word and PDF.
function Pandoc(doc)
  local is_docx = quarto.doc.is_format('docx')
  local is_pdf = quarto.doc.is_format('pdf')
  if not (is_docx or is_pdf) then return nil end
  local authors = doc.meta['by-author']
  local affiliations = doc.meta.affiliations
  if not authors or not affiliations then return nil end

  local names = pandoc.Inlines({})
  for i, author in ipairs(authors) do
    if i > 1 then names:insert(pandoc.Str(',')); names:insert(pandoc.Space()) end
    names:extend(author.name.literal)
    local numbers = {}
    for _, affiliation in ipairs(author.affiliations or {}) do
      numbers[#numbers + 1] = pandoc.utils.stringify(affiliation.number)
    end
    table.sort(numbers, function(a, b) return tonumber(a) < tonumber(b) end)
    if #numbers > 0 then
      names:insert(pandoc.Superscript(pandoc.Inlines({pandoc.Str(table.concat(numbers, ','))})))
    end
  end
  local front = pandoc.Blocks({pandoc.Para(names)})
  for _, affiliation in ipairs(affiliations) do
    local text = pandoc.Inlines({
      pandoc.Superscript(pandoc.Inlines({pandoc.Str(pandoc.utils.stringify(affiliation.number))})),
      pandoc.Space()
    })
    text:extend(affiliation.name)
    front:insert(pandoc.Para(text))
  end
  -- Suppress native bylines only; the source and by-author metadata are retained.
  if is_docx then doc.meta.author = nil end
  if is_pdf then doc.meta.authors = nil end
  front:extend(doc.blocks)
  doc.blocks = front
  return doc
end
