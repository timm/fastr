-- contrast: ranges that most separate two sets of rows.
-- Numeric ranges come from the.bins divisions of -3sd..+3sd.
-- (c) 2026 Tim Menzies <timm@ieee.org> MIT license
local _ENV = setmetatable({}, {__index=require"tbl"})

function discretize(best,rest,out,      f) -- spans, both ways
  f = isNum(best) and _discretizeNum or _discretizeSym
  f(best, rest, out)
  f(rest, best, out) end

function _span(lo,hi,b,r,out) -- discretize's scorer: b^2/(b+r)
  out{lo=lo, hi=hi, score=b*b / (b + r + 1/BIG)} end

function _discretizeSym(best,rest,out,      n0,n1,keys)
  -- one span per symbol
  n0,n1,keys = 1/BIG, 1/BIG, {}
  for k,n in pairs(best) do n0 = n0 + n; keys[k] = true end
  for k,n in pairs(rest) do n1 = n1 + n; keys[k] = true end
  for k in pairs(keys) do
    _span(k, k, (best[k] or 0)/n0, (rest[k] or 0)/n1, out) end
end

function _discretizeNum(best,rest,out, -- the.bins cut spans
                        n,lo,hi,cuts,B,R)
  n  = the.bins
  lo = math.min(best[2]-3*sd(best), rest[2]-3*sd(rest))
  hi = math.max(best[2]+3*sd(best), rest[2]+3*sd(rest))
  cuts,B,R = {},{},{}
  for k = 0,n do
    cuts[k] = lo + k*(hi - lo)/n
    B[k] = norm(best, cuts[k])
    R[k] = norm(rest, cuts[k]) end
  cuts[0], cuts[n] = -BIG, BIG -- extreme spans are open-ended
  for i = 0,n-1 do
    for j = i+1,n do
      if not (i == 0 and j == n) then -- whole column useless
        _span(cuts[i],cuts[j], B[j]-B[i], R[j]-R[i], out) end
    end end end

function contrasts(tbl0,tbl1,      score,all,one,z)
  -- best span over all x columns
  score = function(z1) return z1.score end
  all = most(score)
  for at,_ in pairs(tbl0.cols.x) do
    one = most(score)
    discretize(tbl0.cols.all[at], tbl1.cols.all[at], one)
    z = one()
    if z then
      z.at, z.txt = at, cutTxt(tbl0, at, z)
      all(z) end end
  return all() end

function cutTxt(tbl,at,z,      name) -- span --> readable text
  name = tbl.cols.names[at]
  if not isNum(tbl.cols.all[at]) then
    return ("%s == %s"):format(name, z.lo) end
  if z.lo == -BIG then
    return ("%s <= %.3g"):format(name, z.hi) end
  if z.hi == BIG then
    return ("%s >= %.3g"):format(name, z.lo) end
  return ("%s in %.3g..%.3g"):format(name, z.lo, z.hi) end

function selects(cut,row,      v) -- row inside the range?
  v = row[cut.at]
  return v ~= "?" and cut.lo <= v and v <= cut.hi end


-- --------------------------------------------------------------
function test_discretize() -- separated Nums; span covers one
  local a,b = Num(), Num()
  for _ = 1,100 do add(a, rand()); add(b, rand() + 2) end
  local out = most(function(z) return z.score end)
  discretize(a, b, out)
  local z = out()
  assert(z.score > 0.5 and (z.hi < 2 or z.lo > 1))
  print(("a=N(.5,.3) b=N(2.5,.3): %.2f..%.2f score %.2f")
        :format(math.max(z.lo,-9), math.min(z.hi,9), z.score))
end

function test_contrasts() -- name fastmap's 50:50 top split
  local u = require"unsuper"
  local t = Tbl(csv(the.file))
  local x = function(r1,r2) return u.distx(t,r1,r2) end
  local _,rows = u.fastmap(t.rows, x)
  local half, r0, r1 = #rows // 2, {}, {}
  for i,r in ipairs(rows) do
    if i <= half then r0[#r0+1] = r else r1[#r1+1] = r end end
  local t0, t1 = clone(t,r0), clone(t,r1)
  local z = contrasts(t0, t1)
  local p0, p1 = 0, 0
  for _,r in ipairs(t0.rows) do
    p0 = p0 + (selects(z,r) and 1 or 0) end
  for _,r in ipairs(t1.rows) do
    p1 = p1 + (selects(z,r) and 1 or 0) end
  p0, p1 = p0/#t0.rows, p1/#t1.rows
  print(("best: %s; selects %.0f%% vs %.0f%%")
        :format(z.txt, 100*p0, 100*p1))
  assert(math.abs(p0 - p1) > 0.2) end

function test_halve() -- one line: best split of two halves
  local u = require"unsuper"
  local t = Tbl(csv(the.file))
  local x = function(r1,r2) return u.distx(t,r1,r2) end
  local _,rows = u.fastmap(t.rows, x)
  local half, r0, r1 = #rows // 2, {}, {}
  for i,r in ipairs(rows) do
    if i <= half then r0[#r0+1] = r else r1[#r1+1] = r end end
  local z = contrasts(clone(t,r0), clone(t,r1))
  print(("%-45s %.2f %s"):format(the.file:match"([^/]+)$",
        z and z.score or 0, z and z.txt or "no split")) end

if not pcall(debug.getlocal, 4, 1) then main(_ENV) end
return _ENV
