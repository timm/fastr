-- super: supervised recursive splitting for classify, regress.
-- Each split minimizes the spread of the goals in two halves.
-- (c) 2026 Tim Menzies <timm@ieee.org> MIT license
local _ENV = setmetatable({}, {__index=require"unsuper"})

function disty(tbl,row,      d,n,p,v) -- goals to best corner
  d,n,p = 0, 1/BIG, the.p
  for at,w in pairs(tbl.cols.y) do
    v = norm(tbl.cols.all[at], row[at])
    if v ~= "?" then n,d = n+1, d + math.abs(v - w)^p end end
  return (d/n)^(1/p) end

function yinfo(tbl,      at) -- goal fun, summary fun, spread
  if next(tbl.cols.y) then
    return function(r) return disty(tbl,r) end, Num, sd end
  at = tbl.cols.klass
  return function(r) return r[at] end, Sym, ent end

function split(tbl,rows, -- col+cut minimizing goal spread
               Y,New,div,out,have,ys,left,right,score)
  Y,New,div = yinfo(tbl)
  for at,_ in pairs(tbl.cols.x) do
    if isNum(tbl.cols.all[at]) then
      have = {}
      for _,r in ipairs(rows) do
        if r[at] ~= "?" then have[#have+1] = r end end
      if #have >= 4 then
        have = keysort(have, function(r) return r[at] end)
        ys = {}
        for i,r in ipairs(have) do ys[i] = Y(r) end
        left, right = New(), adds(ys, New())
        for i = 1, #have - 1 do
          add(left, ys[i]); add(right, ys[i], -1)
          if have[i][at] ~= have[i+1][at] then
            score = (i*div(left) +
                     (#have - i)*div(right)) / #have
            if out == nil or score < out.score then
              out = {at=at, score=score,
                     cut=(have[i][at] + have[i+1][at]) / 2}
              out.txt = ("%s <= %.3g"):format(
                          tbl.cols.names[at], out.cut)
            end end end end end end
  return out end

function tree(tbl,rows,      Y,New,go) -- supervised splits
  Y,New = yinfo(tbl)
  function go(rows1,      ys,t,cut,yes,no)
    ys = New()
    for _,r in ipairs(rows1) do add(ys, Y(r)) end
    t = {rows=rows1, n=#rows1, ys=ys}
    if #rows1 > the.stop then
      cut = split(tbl, rows1)
      if cut then
        yes,no = {},{}
        for _,r in ipairs(rows1) do
          if r[cut.at] ~= "?" and r[cut.at] <= cut.cut
          then yes[#yes+1] = r else no[#no+1] = r end end
        if #yes > 0 and #yes < #rows1 then
          t.cut, t.yes, t.no = cut, go(yes), go(no) end end end
    return t end
  return go(rows or tbl.rows) end

function leaf(t,row,      c,v) -- walk row down to its leaf
  while t.cut do
    c,v = t.cut, row[t.cut.at]
    t = (v ~= "?" and v <= c.cut) and t.yes or t.no end
  return t end

function predict(tbl,t,row,k,      rows,at,got) -- knn in leaf
  k = k or the.k
  rows = keysort(leaf(t,row).rows,
           function(r) return distx(tbl,row,r) end)
  at = tbl.cols.klass
  got = Col(tbl.cols.names[at])
  for i = 1, math.min(k,#rows) do add(got, rows[i][at]) end
  return mid(got) end


-- --------------------------------------------------------------
function test_disty() -- sort rows by goals; show best, worst
  local t = Tbl(csv(the.file))
  local lst = keysort(t.rows,
                function(r) return disty(t,r) end)
  assert(disty(t,lst[1]) < disty(t,lst[#lst]))
  print(say(t.cols.names))
  for _,j in ipairs{1,2,3,#lst-2,#lst-1,#lst} do
    print(say(lst[j]), ("d %.2f"):format(disty(t,lst[j]))) end
end

function test_split() -- show split minimizing goal spread
  local t = Tbl(csv(the.file))
  local cut = split(t, t.rows)
  print(cut.txt, ("spread %.3f"):format(cut.score)) end

function test_tree() -- print tree with n and mean goal dist
  local t = Tbl(csv(the.file))
  local t1 = tree(t)
  local n = 0
  for _,l in ipairs(leaves(t1)) do n = n + l.n end
  assert(n == t1.n)
  show(t1, nil, nil,
       function(x) return (" mu=%.2f"):format(mid(x.ys)) end)
  print(#leaves(t1).." leaves") end

function test_classify() -- resubstitution accuracy, diabetes
  local t =
    Tbl(csv"/Users/timm/gits/moot/classify/diabetes.csv")
  local t1, at, acc = tree(t), t.cols.klass, 0
  for _,r in ipairs(t.rows) do
    if predict(t,t1,r) == r[at] then acc = acc + 1 end end
  acc = acc / #t.rows
  print(("resub acc %.2f"):format(acc))
  assert(acc > 0.7) end

if arg[0]:find("super.lua",1,true) then main(_ENV) end
return _ENV
