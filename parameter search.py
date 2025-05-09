from math import *

def refined_search(func, start, end, deltas, lower_limits, upper_limits): # searches for minimum to func, a function taking n integer arguments
  """
  Inputs:
    func - function taking n integers returning a single float
    start - list of n integers (the negative corner of initial search space)
    end - list of n integers (the positive corner of initial search space)
    deltas - list of n integers (step sizes within each argument)
    lower_limits - list of n integers (minimum allowable value for each argument)
    upper_limits - list of n integers (maximum allowable value for each argument)
  Behavior:
    Searches and returns a list of n integers to minimize the value of func.
    The search is done by initially searching the box specified by start and end
    with step sizes given by deltas. Then it repeatedly recenters the box around
    the minimum found, halves both the deltas and the side lengths of the box, and
    searches within the new box. This approach lets it search over a very wide
    sample space in a reasonable amount of time.
  """
  # make copies of start, end, deltas that we can modify
  start = start.copy()
  end = end.copy()
  deltas = deltas.copy()
  
  n = len(start) # length of argument list
  keep_going = True
  while keep_going: # outer loop that manages the boxes
    cur = start.copy() # current parameter list
    best = float('inf') # smallest value of func encountered so far
    best_point = None # parameters that achieved smallest value of func
    while True: # inner loop that searches within the box
      # compute func on cur and modify best and best_point if value is lower
      val = func(*cur)
      if val < best:
        best = val
        best_point = cur.copy()
      # increment cur (in reverse lexicographic order)
      i = 0
      while i < n and cur[i] + deltas[i] > end[i]:
        cur[i] = start[i]
        i += 1
      if i == n:
        break
      cur[i] += deltas[i]
    max_steps = max((end[i] - start[i]) // deltas[i] for i in range(n)) # maximum number of times any coordinate got incremented
    keep_going = False # we will break once the box contains only a single element.
    for i in range(n):
      current_dist = end[i] - start[i] # distance from start[i] to end[i]
      if current_dist == 0: continue
      if deltas[i] > 1: # if deltas[i] > 1, half both deltas[i] and the distance
        deltas[i] //= 2
        new_dist = current_dist // 2
      elif current_dist >= max_steps // 2: # if the distance is roughly the same size as the maximum number of steps, half the distance
        new_dist = current_dist // 2
      else: # if the distance is significantly smaller than the maximum number of steps, leave it the same. This ensures the box won't get too thin.
        new_dist = current_dist
      # recenter the box about best_point
      center = best_point[i]
      half_dist = ceil(new_dist / 2)
      if center-half_dist < lower_limits[i]: # if center is too close to lower_limits[i], start at lower_limits[i].
        start[i] = lower_limits[i]
        end[i] = lower_limits[i] + new_dist
      elif center+half_dist > upper_limits[i]: # if center is too close to upper_limits[i], end at upper_limits[i].
        start[i] = upper_limits[i]-new_dist
        end[i] = upper_limits[i]
      else: # otherwise, start at center-half_dist and end at center-half_dist+new_dist
        start[i] = center - half_dist
        end[i] = start[i] + new_dist # center <= center-half_dist+new_dist <= center+half_dist
      if new_dist > 0: keep_going = True
    print("best:", best, best_point)
  return start

def binary_search(func, start, end, delta=0.00001, sign_preference=True): # binary searches for zero of func, a float-valued function of a single float argument
  """
  Inputs:
    func - function taking float and returning float
    start - float (lower endpoint of search interval)
    end - float (upper endpoint of search interval)
    delta - float (desired absolute error)
    sign_preference - bool (the desired value of func(x) > 0, where x is the return value)
  Behavior:
    Binary searches to find zero of func within interval [start, end].
    Returns x within delta of the zero such that (func(x) > 0) == sign_preference
  """
  stval = func(start)
  enval = func(end)
  if stval == 0:
    return start
  if enval == 0:
    return end
  ssign = stval > 0
  esign = enval > 0
  assert start < end and ssign != esign # assert that start < end and func takes different signs on start and end
  while end-start > delta:
    mid = (start + end) / 2
    msign = func(mid) > 0
    if msign == ssign:
      start = mid
    else:
      end = mid
  if ssign == sign_preference:
    return start
  else:
    return end

tlst1 = 2 * log(sqrt(2) + 1) # frequently used constant

def eq_5_32_value(Kn, L, R1, R2, mun, verbose=False): # system for y = infinity
  """
  p is a global variable representing the p in the equation x^2-2 = y^p
  Inputs:
    Kn - integer, that is 100 times the value of K', so that we can search on the integer values
    L - integer
    R1 - integer
    R2 - integer
    mun - integer, that is 100 times the value of mu, so that we can search on the integer values
  Behavior:
    Returns smallest value of the quantity log(rho) mu L K' from equation (5.32)
    with the given choices of K', L, R1, R2, mu
  """
  if Kn < 1 or L < 1 or R1 < 1 or R2 < 1:
    return float('inf')
  # set K' and mu from Kn and mun
  Kp = Kn / 100 # K'
  mu = mun / 100
  # compute limiting values of R, C2, log(b), g, sigma
  R = R1 + R2 - 1
  C2 = Kp * L / min(R2, p)
  log_b = log(L * p / (2 * min(R2, p))) + 3/2
  g = 1/4 - min(R2, p) / (12 * R)
  sigma = (1 + 2 * mu - mu**2) / 2
  # compute A, B, C so that constraint (5.18) is A log(rho) - B rho > C
  A = Kp * (sigma * L - 1)
  B = g * L * tlst1 * C2
  C = 2 * Kp * log_b + g * L * (2 * R + tlst1 * C2)

  if verbose:
    print(f"A, B, C: {A}, {B}, {C}")
  func = lambda rho: A * log(rho) - B * rho - C # find smallest rho >= 1 that satisfies func(rho) > 0
  start = 1 # lower limit for rho
  end = A / B # maximum of func
  if start < end and func(start) <= 0 and func(end) > 0: # rho is between start and end
    rho = binary_search(func, start, end)
  elif func(start) > 0: # rho = 1 works
    rho = 1
  else: # no value of rho works (either end <= 1 and func(start) <= 0 or func(end) <= 0)
    return float('inf')
  if verbose:
    print(f'rho: {rho}')
  return mu * log(rho) * Kp * L

start = [3000, 5, 1, 100, 33]
end = [8000, 15, 10, 150, 100]
deltas = [100, 2, 2, 5, 5]
lower_limits = [1, 1, 1, 1, 33]
upper_limits = [1000000, 10000, 10000, 10000, 100]

def test_p(p_):
  global p
  p = p_
  Kn, L, R1, R2, mun = refined_search(eq_5_32_value, start, end, deltas, lower_limits, upper_limits)
  return eq_5_32_value(Kn, L, R1, R2, mun) < p / 2

print(test_p(1832))
