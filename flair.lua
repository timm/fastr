local lib = require"lib"
local _ENV = setmetatable({}, {__index=lib})
the,help = {},[[
start: config and a tiny test-runner CLI.
(c) 2026 Tim Menzies <timm@ieee.org> MIT license

Options:

  -h             show help
  --p=2          minkowski coefficient
  --stop=4       stopping rule for recursive tree generation
  --few=256      sub-sample size for pole finding
  --k=5          nearest neighbors used in a leaf
  --check=5      optimization: how many top picks to evaluate
  --bins=10      number of bins for discretization
  --cliffs=.197  cliffs delta: max effect size for "same"
  --cohen=.35    cohen d: max mean separation for "same"
  --seed=1234    random number generation
  --file=/Users/timm/gits/moot/optimize/misc/auto93.csv]]

function atom(s) -- string --> number or trimmed string
  if type(s) ~= "string" then return s end
  return tonumber(s) or s:match"^%s*(.-)%s*$" end

for k,v in help:gmatch("(%w+)=(%S+)") do the[k]=atom(v) end

function csv(file,      f) -- iterate a csv file's atom rows
  f = assert(io.open(file))
  return function(    s,t) 
    s = f:read()
    if s == nil then f:close() else
      t={}; for x in s:gmatch"([^,]+)" do t[#t+1]=atom(x) end
      return t end end end

Cols,Num,Sym = {},{},{}

function Sym.new() return {} end
function Num.new() return {n=0, mu=0, m2=0} end
function Col(s)    return (s:find"^%u" and Num or Sym).new() end
function Tbl(src)  return adds(src, {rows={}, cols=nil}) end

function Cols.new(names,     i,__roles)
  i = {names=names, all={}, x={}, y={}, klass=nil}
  for at,name in ipairs(names) do
    i.all[at] = Col(name)
    Cols.roles(name:sub(-1), at) end
  return i end

function Cols.roles(z,at)
  if     z == "X" then return
  elseif z == "!" then i.klass = at
  elseif z == "+" then i.y[at] = 1
  elseif z == "-" then i.y[at] = 0
  else   i.x[at] = at end end 

function adds(src,i) -- add all from a list or iterator
  i = i or Num()
  if type(src) == "table"
  then for _,x in ipairs(src) do add(i,x) end
  else for x in src           do add(i,x) end end
  return i end

function add(i,v,  w)
  w = w or 1
  if     i.mu then welford(i,v,w)
  elseif i.all then for j,c in pairs(i.all) do add(c,v[j],w) end
  elseif i.rows then 
    if   i.cols 
    then i.rows[1+#i.rows]=v; add(i.cols, v,   w)
    else i.cols = Cols.new(v) end
  else   i[v] = w + (i[v] or 0) end
  return v end

function welford(i,v,w)
  i.n  = i.n + w
  if i.n <= 1 then i.n, i.mu, i.m2 = 0,0,0 else 
    d    = v - i.mu
    i.mu = i.mu + w*d/i.n
    i.m2 = i.m2 + w*d*(v - i.mu) end end

function sort(t,f) table.sort(t,f); return t end

floor,cat = math.floor,table.concat

function say(t,      u) -- anything --> string, tidy numbers
  if type(t)=="number" then
    return (t==floor(t) and "%d" or "%.3f"):format(t) end
  if type(t)~="table" then return tostring(t) end
  u={}
  for k,v in pairs(t) do
    u[#u+1] = (#t>0 and "" or k.."=")..say(v) end
  return "{"..cat(#t==0 and sort(u) or u,", ").."}" end

SEED = the.seed
function rand(lo,hi) -- pseudo-random lo..hi (default 0..1)
  lo, hi = lo or 0, hi or 1
  SEED = (16807 * SEED) % 2147483647
  return lo + (hi - lo) * SEED / 2147483647 end

function rint(lo,hi) -- pseudo-random integer lo..hi
  return math.floor(0.5 + rand(lo,hi)) end

function run(f,      ok,err) -- reseed, call f, catch crashes
  SEED = the.seed
  ok,err = pcall(f)
  if not ok then print(err) end end

eg = eg or {} -- demo table: eg["-x"] = function(v) ... end

if arg[0] and arg[0]:find"flair" then
  for j,s in ipairs(arg) do
    if eg[s] then eg[s](arg[j+1]) end
    s = s:gsub("^[-]+","")
    if the[s] ~= nil then the[s]=atom(arg[j+1]) end end end 

function demos(env,      t,k,doc) -- local test_ funs, in order
  t = {}
  for line in io.lines(arg[0]) do
    k,doc = line:match"^function%s+(test_[%w_]+).-[-][-]%s*(.*)"
    if k and rawget(env,k) then t[#t+1] = {k,doc} end end
  return t end

function main(env,      k,eq,v,f) -- -demo | --setting value
  for i,s in ipairs(arg) do
    k,eq,v = s:match"^[-][-](%w+)(=?)(.*)"
    if s == "-h" then
      print(help.."\n\nDemos:\n")
      for _,kf in ipairs(demos(env)) do
        print(("  -%-13s %s"):format(kf[1]:sub(6), kf[2])) end
    elseif s == "-all" then
      for _,kf in ipairs(demos(env)) do
        print("\n# "..kf[1]); run(env[kf[1]]) end
    elseif k and the[k] ~= nil then
      the[k] = atom(eq == "=" and v or arg[i+1])
    elseif s:find"^[-]" then
      f = rawget(env, "test_"..s:sub(2):gsub("-","_"))
      if f then run(f) end end end end


-- --------------------------------------------------------------
function test_the() -- show current settings
  print(say(the)) end

function test_atom() -- strings coerce to numbers or strings
  assert(atom"2" == 2 and atom"2.1" == 2.1)
  assert(atom" a " == "a")
  print("'2' ->", atom"2", "| '2.1' ->", atom"2.1",
        "| ' a ' ->", "'"..atom" a ".."'") end

function test_csv() -- csv reader finds many rows in the.file
  local n = 0
  for _ in csv(the.file) do n = n + 1 end
  assert(n > 100)
  print(n, "rows") end

if arg[0]:find("start.lua",1,true) then main(_ENV) end
return _ENV
