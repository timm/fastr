-- unsuper: unsupervised recursive bi-clustering (no goals).
-- Fastmap poles halve the data; contrast names the split.
-- (c) 2026 Tim Menzies <timm@ieee.org> MIT license
local _ENV = setmetatable({}, {__index=require"contrast"})

function distx(tbl,r1,r2,      d,n,p) -- x-column distance
  d,n,p = 0, 1/BIG, the.p
  for at,_ in pairs(tbl.cols.x) do
    n = n + 1
    d = d + _distx(tbl.cols.all[at], r1[at], r2[at])^p end
  return (d/n)^(1/p) end

function _distx(col,a,b) -- helper for one column
  if a == "?" and b == "?" then return 1 end
  if not isNum(col) then return a == b and 0 or 1 end
  a, b = norm(col,a), norm(col,b)
  if a == "?" then a = b < 0.5 and 1 or 0 end
  if b == "?" then b = a < 0.5 and 1 or 0 end
  return math.abs(a - b) end

function proj(x,lohi,row,      a,b) -- project row on a line
  a,b = x(row,lohi.lo), x(row,lohi.hi)
  return (a*a + lohi.c*lohi.c - b*b) / (2*lohi.c + 1/BIG) end

function fastmap(rows,x,      some,w,a,b,m1,m2,lohi)
  -- sort rows on a line between 2 distant poles
  some = many(rows, math.min(the.few, #rows))
  w = any(some)
  m1 = most(function(r) return x(w,r) end)
  for _,r in ipairs(some) do m1(r) end
  a = m1()
  m2 = most(function(r) return x(a,r) end)
  for _,r in ipairs(some) do m2(r) end
  b = m2()
  lohi = {lo=a, hi=b, c=x(a,b)}
  return lohi, keysort(rows,
           function(r) return proj(x,lohi,r) end) end

function cluster(tbl,rows,      x,go)
  -- bi-cluster; contrast names the splits
  x = function(r1,r2) return distx(tbl,r1,r2) end
  function go(rows1,      t,tmp,half,r0,r1,cut,yes,no)
    t = {rows=rows1, n=#rows1}
    if #rows1 > the.stop then
      tmp = select(2, fastmap(rows1,x))
      half, r0, r1 = #tmp // 2, {}, {}
      for i,r in ipairs(tmp) do
        if i <= half then r0[#r0+1]=r else r1[#r1+1]=r end end
      cut = contrasts(clone(tbl,r0), clone(tbl,r1))
      if cut then
        yes,no = {},{}
        for _,r in ipairs(rows1) do
          if selects(cut,r) then yes[#yes+1] = r
          else no[#no+1] = r end end
        if #yes > 0 and #yes < #rows1 then
          t.cut, t.yes, t.no = cut, go(yes), go(no) end end end
    return t end
  return go(rows or tbl.rows) end

function show(t,pre,txt,more,      sub) -- print tree with n
  more = more or function(_) return "" end
  print((pre or "")..(txt or "").."n="..t.n..more(t))
  if t.cut then
    sub = pre == nil and "" or pre.."|  "
    show(t.yes, sub, t.cut.txt.."; ", more)
    show(t.no,  sub, "!"..t.cut.txt.."; ", more) end end

function leaves(t,out) -- all terminal nodes of a tree
  out = out or {}
  if t.cut then leaves(t.yes,out); leaves(t.no,out)
  else out[#out+1] = t end
  return out end


-- --------------------------------------------------------------
function test_distx() -- symmetric, zero on self, in 0..1
  local t = Tbl(csv(the.file))
  local r1, r2 = t.rows[1], t.rows[#t.rows]
  assert(distx(t,r1,r1) == 0)
  assert(distx(t,r1,r2) == distx(t,r2,r1))
  for _,r in ipairs(t.rows) do
    local d = distx(t,r1,r)
    assert(d >= 0 and d <= 1) end
  print(("d(r1,r1) %.2f | d(r1,r2) %.2f")
        :format(distx(t,r1,r1), distx(t,r1,r2)))
  print("r1 "..cat(r1).."\nr2 "..cat(r2)) end

function test_fastmap() -- rows sorted along a line of 2 poles
  local t = Tbl(csv(the.file))
  local x = function(r1,r2) return distx(t,r1,r2) end
  local lohi, rows = fastmap(t.rows, x)
  assert(lohi.c > 0)
  local last = -BIG
  for _,r in ipairs(rows) do
    local p = proj(x,lohi,r)
    assert(p >= last); last = p end
  print("pole1 "..cat(lohi.lo).."\npole2 "..cat(lohi.hi))
  print(("separation %.2f"):format(lohi.c)) end

function test_cluster() -- print tree; leaf counts sum to n
  local t = Tbl(csv(the.file))
  local tr = cluster(t)
  assert(tr.n == #t.rows)
  local n = 0
  for _,l in ipairs(leaves(tr)) do n = n + l.n end
  assert(n == tr.n)
  show(tr)
  print(#leaves(tr).." leaves") end

if not pcall(debug.getlocal, 4, 1) then main(_ENV) end
return _ENV
