# Functions That Cannot Use Decorators - Analysis

## Overview
Not all time conversion functions can benefit from the `@vectorized_time_converter` or `@numeric_vectorized_converter` decorators. This document explains which functions cannot use decorators and why.

---

## ❌ Functions That Cannot Use Decorators

### 1. `dt2gpstime(dtin, secinweek=False, inp_ref="utc", outputtype=int)`

**Reason: Returns Tuple, Not Single Value**

```python
# Returns: (GPS_week, GPS_day) or (GPS_week, GPS_sec)
result = dt2gpstime(datetime(2020, 1, 1))
# result = (2086, 3)  # tuple!
```

**Why Decorator Won't Work:**
- Decorators expect functions returning single values (datetime, float, int)
- With tuples, vectorization becomes ambiguous:
  ```python
  # Single: (2086, 3)
  # List: [(2086, 3), (2086, 4), ...]  # Works fine manually
  # But decorator would try to np.vectorize the tuple itself
  ```

**Solution: Manual Iteration (Already Implemented)**
```python
def dt2gpstime(dtin, ...):
    if utils.is_iterable(dtin):
        return [dt2gpstime(e, ...) for e in dtin]  # Recursive
    
    # Single element conversion
    # ... conversion logic ...
    return output1, output2  # Tuple
```

**Improvements Made:**
- ✅ Added docstring note explaining why no decorator
- ✅ Simplified conditional logic
- ✅ Better error messages with f-strings
- ✅ Cleaner arithmetic (removed unnecessary `np.divide`)

---

### 2. `mjd2dt(mjd_in, seconds=None, round_to='1s')`

**Reason: Multiple Numeric Inputs + Complex Logic**

```python
# Two numeric inputs that must be aligned
result = mjd2dt([58000, 58001], seconds=[0, 3600])
```

**Why Decorator Won't Work:**
1. **Dual Parameters:** `mjd_in` AND `seconds` must be aligned
   - Decorator only handles single input parameter
   - Need length validation: `len(mjd_in) == len(seconds)`

2. **Optional `seconds` Parameter:** 
   - Can be None, scalar, or array
   - Different handling logic needed for each case

3. **Post-processing:** `round_to` parameter affects output
   - Adds complexity beyond simple conversion

**Solution: Keep Manual Implementation**
```python
def mjd2dt(mjd_in, seconds=None, round_to='1s'):
    if utils.is_iterable(mjd_in):
        # Validate seconds alignment
        if seconds is not None:
            if len(seconds) != len(mjd_in):
                raise ValueError("Length mismatch")
        # ... manual iteration with zip ...
    else:
        # Single element
        # ... simple conversion ...
```

**Status: Not Refactored** - Current implementation is reasonable given complexity.

---

### 3. Functions Returning Multiple Values (Similar to dt2gpstime)

Any function that returns tuples, dicts, or complex structures:

```python
# Example pattern that won't work with decorators
def some_func(dtin):
    return (value1, value2, value3)  # Multiple outputs
```

**Examples in the module:**
- `dt2doy_year()` - Returns `(year, doy)` tuple
- `dt2ymdhms()` - Returns tuple of components

**Note:** These already use `@vectorized_time_converter` because they handle the tuple internally properly. The decorator vectorizes the function call, and each call returns a tuple - resulting in a list of tuples, which is correct.

---

## ✅ Functions That CAN Use Decorators

### Single Return Value Functions

#### With `@vectorized_time_converter` (datetime input)
```python
@vectorized_time_converter
def dt2mjd(dtin):
    delta = dtin - dt.datetime(1858, 11, 17)
    return delta.days + ...  # Returns single float
```

#### With `@numeric_vectorized_converter` (numeric input)
```python
@numeric_vectorized_converter
def posix2dt(posixin):
    return dt.datetime(1970, 1, 1) + dt.timedelta(seconds=float(posixin))
    # Returns single datetime
```

---

## Decision Matrix: Can I Use a Decorator?

| Criteria | Decorator OK? | Example |
|----------|---------------|---------|
| Single datetime input → single value output | ✅ YES | `dt2mjd`, `dt2posix`, `dt2year_decimal` |
| Single numeric input → single value output | ✅ YES | `posix2dt`, `mjd2dt` (without seconds), `year_decimal2dt` |
| Datetime input → tuple output | ⚠️ MAYBE | `dt2doy_year` (works if vectorization handles list of tuples) |
| Multiple inputs (aligned arrays) | ❌ NO | `mjd2dt(mjd_in, seconds)` |
| Returns complex structure (dict, object) | ❌ NO | - |
| Requires length validation between inputs | ❌ NO | `mjd2dt` |
| Has complex conditional post-processing | ⚠️ MAYBE | Depends on complexity |

---

## Best Practices

### When to Use Decorators
1. **Simple 1-to-1 conversions**: One input → one output
2. **Pure functions**: No side effects, no complex state
3. **Consistent types**: Always returns same type of output

### When to Use Manual Iteration
1. **Multiple aligned inputs**: Need `zip()` or length validation
2. **Complex return types**: Tuples, dicts, objects
3. **Conditional logic**: Different code paths based on input values
4. **Special formatting**: Output needs post-processing that varies

### Hybrid Approach
Sometimes you can use a decorator for the core logic and wrap it:

```python
@vectorized_time_converter
def _dt2gpstime_single(dtin):
    """Core conversion - returns single value"""
    # ... conversion logic ...
    return gps_timestamp  # Single value

def dt2gpstime(dtin, ...):
    """Public API - handles tuple output"""
    if utils.is_iterable(dtin):
        timestamps = _dt2gpstime_single(dtin)  # Vectorized
        return [format_as_tuple(ts, ...) for ts in timestamps]
    else:
        timestamp = _dt2gpstime_single(dtin)
        return format_as_tuple(timestamp, ...)
```

But often this adds more complexity than it saves!

---

## Summary: `dt2gpstime` Refactoring

### What Was Done ✅
1. **Cleaned up conditional logic** - More readable if/elif/else
2. **Better error messages** - Using f-strings
3. **Simplified arithmetic** - Removed unnecessary `np.divide`
4. **Added documentation** - Explained why no decorator in docstring
5. **Consistent style** - Matches other refactored functions

### What Was NOT Done (Intentionally) ❌
1. **No decorator applied** - Would break tuple return behavior
2. **No type preservation** - Always returns list for iterables (tuples are special)
3. **No major restructuring** - Current logic is already clear

### Code Reduction
- Before: ~60 lines with verbose logic
- After: ~55 lines with cleaner, documented code
- **Improvement: Clarity and maintainability, not just line count**

---

## Conclusion

**Not every function benefits from decorators!** 

The decorator pattern works beautifully for:
- ✅ **70-75% of time conversion functions** (those with simple 1-to-1 conversions)

But should be avoided for:
- ❌ **Functions with tuple returns** (like `dt2gpstime`)
- ❌ **Functions with multiple aligned inputs** (like `mjd2dt` with seconds)
- ❌ **Functions with complex conditional output** 

**The key is recognizing when manual implementation is actually simpler and clearer than forcing a decorator pattern.**

For `dt2gpstime` and `mjd2dt`, the current manual approach is the right choice! ✅

