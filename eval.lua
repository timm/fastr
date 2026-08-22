-- eval: cross-val, confusion, effect size, significance, ranks.
-- (c) 2026 Tim Menzies <timm@ieee.org> MIT license
local _ENV = setmetatable({}, {__index=require"super"})

function confuse(pairs_,      out,got,want) -- per-class stats
  out = {}
  for _,p in ipairs(pairs_) do
    for i = 1,2 do
      out[p[i]] = out[p[i]] or
        {label=p[i], tp=0, fp=0, fn=0, tn=0} end end
  for _,p in ipairs(pairs_) do
    got,want = p[1], p[2]
    for k,c in pairs(out) do
      if k == want then
        if got == want then c.tp=c.tp+1 else c.fn=c.fn+1 end
      elseif k == got then c.fp = c.fp + 1
      else c.tn = c.tn + 1 end end end
  for _,c in pairs(out) do
    c.acc  = (c.tp+c.tn) / (c.tp+c.fn+c.fp+c.tn + 1/BIG)
    c.pd   = c.tp / (c.tp + c.fn + 1/BIG)
    c.pf   = c.fp / (c.fp + c.tn + 1/BIG)
    c.prec = c.tp / (c.tp + c.fp + 1/BIG) end
  return out end

function xval(rows,n,m,fun,      tmp,train,test) -- n*m folds
  n, m = n or 5, m or 5
  tmp = copy(rows)
  for _ = 1,m do
    shuffle(tmp)
    for j = 1,n do
      train,test = {},{}
      for i,r in ipairs(tmp) do
        if i % n == j - 1 then test[#test+1] = r
        else train[#train+1] = r end end
      fun(train,test) end end end

function cliffs(xs,ys,      gt,lt) -- effect size small?
  gt,lt = 0,0
  for _,x in ipairs(xs) do
    for _,y in ipairs(ys) do
      if x > y then gt = gt + 1 end
      if x < y then lt = lt + 1 end end end
  return math.abs(gt - lt) / (#xs * #ys + 1/BIG)
         <= the.cliffs end

function cohen(xs,ys,      a,b,sp) -- means differ < d*sd?
  a,b = adds(xs), adds(ys)
  sp = (((a[1]-1)*sd(a)^2 + (b[1]-1)*sd(b)^2) /
        (a[1] + b[1] - 2 + 1/BIG))^0.5
  return math.abs(a[2] - b[2]) <= the.cohen * sp end

function ks(xs,ys,      cdf,d) -- kolmogorov-smirnov same?
  function cdf(t,v,      n)
    n = 0
    for _,x in ipairs(t) do if x <= v then n = n + 1 end end
    return n / #t end
  d = 0
  for _,t in ipairs{xs,ys} do
    for _,v in ipairs(t) do
      d = math.max(d, math.abs(cdf(xs,v) - cdf(ys,v))) end end
  return d <= 1.36 * ((#xs + #ys) / (#xs * #ys))^0.5 end

function same(xs,ys) -- statistically indistinguishable?
  return cliffs(xs,ys) and ks(xs,ys) end

function ranks(d,      names,out,rank,last) -- same() = tied
  names = {}
  for k in pairs(d) do names[#names+1] = k end
  names = keysort(names,
            function(k) return adds(d[k])[2] end)
  out, rank = {}, 0
  for _,k in ipairs(names) do
    if last and not same(last, d[k]) then rank = rank + 1 end
    out[k], last = rank, d[k] end
  return out end


-- --------------------------------------------------------------
function test_confuse() -- confusion stats on a small example
  local ps = {}
  for _ = 1,4 do ps[#ps+1] = {"a","a"} end
  ps[#ps+1] = {"b","a"}
  for _ = 1,2 do ps[#ps+1] = {"b","b"} end
  ps[#ps+1] = {"a","b"}
  local a = confuse(ps).a
  assert(a.tp == 4 and a.fn == 1 and a.fp == 1 and a.tn == 2)
  for _,c in pairs(confuse(ps)) do
    print(("%s acc %.2f pd %.2f pf %.2f prec %.2f")
          :format(c.label, c.acc, c.pd, c.pf, c.prec)) end end

function test_stats() -- same() ok on tiny shifts, not big
  local a,b,c = {},{},{}
  for i = 1,40 do
    a[i] = rand(); b[i] = a[i] + 0.05; c[i] = a[i] + 2 end
  assert(same(a,b) and not same(a,c))
  assert(cohen(a,b) and not cohen(a,c))
  print("a=b: "..tostring(same(a,b))..
        "  a=c: "..tostring(same(a,c))) end

function test_ranks() -- same rank iff same distribution
  local d = {x1={}, x2={}, x3={}, x4={}}
  for i = 1,30 do
    d.x1[i] = rand()
    d.x2[i] = rand()
    d.x3[i] = rand() + 0.5
    d.x4[i] = rand() + 2 end
  local r = ranks(d)
  print(say(r))
  assert(r.x1 == r.x2 and r.x2 < r.x3 and r.x3 < r.x4) end

function test_xval() -- cross-val: diabetes and soybean
  for _,f in ipairs{"diabetes","soybean"} do
    print("\n## "..f)
    local t =
      Tbl(csv("/Users/timm/gits/moot/classify/"..f..".csv"))
    local k, ps = t.cols.klass, {}
    xval(t.rows, 3, 5, function(train,test)
      local tr = clone(t, train)
      local t1 = tree(tr)
      for _,row in ipairs(test) do
        ps[#ps+1] = {predict(tr,t1,row), row[k]} end end)
    local cs = {}
    for _,c in pairs(confuse(ps)) do cs[#cs+1] = c end
    for _,c in ipairs(keysort(cs,
                        function(c) return -c.pd end)) do
      print(("%25s acc %.2f pd %.2f pf %.2f prec %.2f")
            :format(c.label, c.acc, c.pd, c.pf, c.prec)) end end
end

function test_opt() -- optimizer: tree picks beat random
  local t = Tbl(csv(the.file))
  local best, rnd = Num(), Num()
  xval(t.rows, 3, 10, function(train,test)
    local tr = clone(t, train)
    local t1 = tree(tr)
    local picks = keysort(test,
      function(r) return mid(leaf(t1,r).ys) end)
    local b, r = BIG, BIG
    for i = 1, the.check do
      b = math.min(b, disty(tr, picks[i]))
      r = math.min(r, disty(tr, any(test))) end
    add(best, b); add(rnd, r) end)
  local lo = BIG
  for _,r in ipairs(t.rows) do
    lo = math.min(lo, disty(t,r)) end
  print(("true best %.2f | best of %s picked %.2f | random %.2f")
        :format(lo, the.check, mid(best), mid(rnd))) end

if arg[0]:find("eval.lua",1,true) then main(_ENV) end
return _ENV
