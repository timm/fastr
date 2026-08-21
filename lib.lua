-- lib: general helpers (random, lists, strings, csv).
-- (c) 2026 Tim Menzies <timm@ieee.org> MIT license
local _ENV = setmetatable({}, {__index=_G})

BIG  = 1e32
SEED = 1234

function atom(s) -- string --> number or trimmed string
  if type(s) ~= "string" then return s end
  return tonumber(s) or s:match"^%s*(.-)%s*$" end

function rand(lo,hi) -- pseudo-random lo..hi (default 0..1)
  lo, hi = lo or 0, hi or 1
  SEED = (16807 * SEED) % 2147483647
  return lo + (hi - lo) * SEED / 2147483647 end

function rint(lo,hi) -- pseudo-random integer lo..hi
  return math.floor(0.5 + rand(lo,hi)) end

function any(t) -- pick one item at random
  return t[rint(1,#t)] end

function many(t,n,      u) -- pick n items at random
  u = {}
  for _ = 1,n do u[#u+1] = any(t) end
  return u end

function shuffle(t,      j) -- randomize order, in place
  for i = #t,2,-1 do
    j = rint(1,i)
    t[i],t[j] = t[j],t[i] end
  return t end

function copy(t,      u) -- shallow copy of a list
  u = {}
  for i,v in ipairs(t) do u[i] = v end
  return u end

function keysort(t,fun,      tmp,u) -- sorted copy, order fun
  tmp = {}
  for i,v in ipairs(t) do tmp[i] = {fun(v), v} end
  table.sort(tmp, function(a,b) return a[1] < b[1] end)
  u = {}
  for i,p in ipairs(tmp) do u[i] = p[2] end
  return u end

function most(fun,      b) -- carried max: f(x) adds, f() gets
  return function(x)
    if x ~= nil and (b == nil or fun(x) > fun(b)) then b = x end
    return b end end

function cat(t,      u,n,ks) -- anything --> string
  if type(t) == "number" then return ("%g"):format(t) end
  if type(t) ~= "table"  then return tostring(t) end
  u, n = {}, 0
  for _ in pairs(t) do n = n + 1 end
  if n == #t then
    for _,v in ipairs(t) do u[#u+1] = cat(v) end
  else
    ks = {}
    for k in pairs(t) do ks[#ks+1] = k end
    table.sort(ks, function(a,b)
      return tostring(a) < tostring(b) end)
    for _,k in ipairs(ks) do
      u[#u+1] = tostring(k).."="..cat(t[k]) end end
  return "{"..table.concat(u, ", ").."}" end

function csv(file,      f) -- iterate a csv file's atom rows
  f = assert(io.open(file))
  return function() -- (for-in iterator: locals must stay local)
    local s = f:read()
    if s == nil then f:close() else
      local t = {}
      for x in s:gmatch"([^,]+)" do t[#t+1] = atom(x) end
      return t end end end

return _ENV
