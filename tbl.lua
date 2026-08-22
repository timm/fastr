-- tbl: Num, Sym, Cols, Tbl: incremental column summaries.
-- (c) 2026 Tim Menzies <timm@ieee.org> MIT license
local _ENV = setmetatable({}, {__index=require"start"})

NUM = {} -- metatable marking Nums

function Num() -- {n, mu, m2}: tracks a column of numbers
  return setmetatable({0,0,0}, NUM) end

function Sym() -- {symbol=count}: tracks a column of symbols
  return {} end

function isNum(col) -- is this column a Num?
  return getmetatable(col) == NUM end

function Col(s) -- Nums start with an upper case letter
  return s:find"^%u" and Num() or Sym() end

function Cols(names,      i,z) -- names --> x,y columns
  i = {names=names, all={}, x={}, y={}, klass=nil}
  for at,name in ipairs(names) do
    i.all[at] = Col(name)
    z = name:sub(-1)
    if z ~= "X" then
      if     z == "!" then i.klass = at
      elseif z == "+" then i.y[at] = 1
      elseif z == "-" then i.y[at] = 0
      else   i.x[at] = true end end end
  return i end

function Tbl(src,      i) -- row iterator (first=names) --> tbl
  for row in src do
    if i then add(i,row)
    else i = {rows={}, cols=Cols(row)} end end
  return i end

function clone(tbl,rows,      i) -- empty copy, same columns
  i = {rows={}, cols=Cols(tbl.cols.names)}
  for _,r in ipairs(rows or {}) do add(i,r) end
  return i end

function add(i,v,w) -- add value (or row) v, weight w
  w = w or 1
  if v ~= "?" then
    if     isNum(i) then welford(i,v,w)
    elseif i.rows   then i.rows[#i.rows+1]=v; add(i.cols,v,w)
    elseif i.all    then
      for at,col in ipairs(i.all) do add(col, v[at], w) end
    else   i[v] = (i[v] or 0) + w end end
  return i end

function adds(src,i,f) -- src, f: anything items() accepts
  i, f = i or Num(), items(f)
  for x in items(src) do add(i, f(x)) end
  return i end

function welford(i,v,w,      n,d) -- update a Num in place
  n = i[1] + w
  if n < 1 then i[1],i[2],i[3] = 0,0,0 else
    d = v - i[2]
    i[2] = i[2] + w*d/n
    i[3] = i[3] + w*d*(v - i[2])
    i[1] = n end end

function mid(col,      hi,mode) -- middle of a distribution
  if isNum(col) then return col[2] end
  hi = -1
  for k,n in pairs(col) do
    if n > hi then hi,mode = n,k end end
  return mode end

function mids(tbl,      u) -- middles of all columns
  u = {}
  for at,col in ipairs(tbl.cols.all) do u[at] = mid(col) end
  return u end

function sd(num) -- diversity of a Num
  return num[1] <= 1 and 0 or (num[3]/(num[1]-1))^0.5 end

function ent(d,      N,e) -- diversity of a Sym
  N,e = 0,0
  for _,n in pairs(d) do N = N + n end
  for _,n in pairs(d) do
    if n > 0 then e = e - n/N * math.log(n/N, 2) end end
  return e end

function norm(col,x,      z) -- x --> 0..1, logistic cdf
  if x == "?" then return x end
  z = (x - col[2]) / (sd(col) + 1/BIG)
  z = math.max(-3, math.min(3, z))
  return 1/(1 + math.exp(-1.7 * z)) end


-- --------------------------------------------------------------
function test_num() -- Num tracks n,mu,sd; subtract restores
  local i = Num()
  for _ = 1,100 do add(i, rand()) end
  local n,mu = i[1],i[2]
  assert(n == 100 and math.abs(mu - 0.5) < 0.1)
  assert(math.abs(sd(i) - 12^-0.5) < 0.05)
  add(add(i, 0.5), 0.5, -1)
  assert(math.abs(i[2] - mu) < 1e-9 and i[1] == n)
  print(("100 rands: n %s mu %.3f sd %.3f"):format(n,mu,sd(i)))
end

function test_sym() -- Sym counts symbols; mid is the mode
  local i = Sym()
  for c in ("aabbbc"):gmatch"." do add(i,c) end
  assert(i.b == 3 and mid(i) == "b")
  print(say(i), "mid: "..mid(i)) end

function test_tbl() -- load the.file; find rows and x,y columns
  local t = Tbl(csv(the.file))
  assert(#t.rows > 100 and next(t.cols.x) and next(t.cols.y))
  print(#t.rows.." rows; x: "..say(t.cols.x)..
        " y: "..say(t.cols.y)) end

function test_mid() -- show middle of every column in the.file
  local t = Tbl(csv(the.file))
  local m = mids(t)
  assert(#m == #t.cols.names)
  print(say(m)) end

function test_norm() -- norm maps to 0..1, monotonically
  local num = Num()
  for _ = 1,100 do add(num, rand(0,1000)) end
  local last = 0
  for x = 0,999,99 do
    local v = norm(num,x)
    assert(v >= last and v >= 0 and v <= 1); last = v end
  local u = {}
  for _,x in ipairs{0,250,500,750,999} do
    u[#u+1] = ("%d->%.2f"):format(x, norm(num,x)) end
  print(table.concat(u," ")) end

function test_cdf() -- norm(hi)-norm(lo) ~ mass inside lo..hi
  local xs, num = {}, Num()
  for i = 1,1000 do
    local x = -6
    for _ = 1,12 do x = x + rand() end
    xs[i] = x; add(num,x) end
  for _,lohi in ipairs{{-1,1},{0,1},{-2,2}} do
    local lo,hi = lohi[1],lohi[2]
    local got = norm(num,hi) - norm(num,lo)
    local want = 0
    for _,x in ipairs(xs) do
      if lo <= x and x <= hi then want = want + 1 end end
    want = want / #xs
    print(("mass %2s..%s: cdf %.2f truth %.2f")
          :format(lo,hi,got,want))
    assert(math.abs(got - want) < 0.05) end end

if arg[0]:find("tbl.lua",1,true) then main(_ENV) end
return _ENV
