# Metasmith Result Compilation Bug Report

## Summary

The KofamScan E2E test fails during result compilation in `metasmith/agents.py:1016` with:
```
AssertionError: assert len(to_del)>0
```

The Nextflow workflow completes successfully and produces valid output, but the metasmith agent crashes when compiling results.

## Root Cause Analysis

### The Bug

In `RunWorkflow()` (agents.py lines 910-1020), during result compilation:

1. **`path2inst`** is built by iterating over staged data libraries using `lib.Iterate()` (line 928-930)
2. **`kv2path`** is built from input manifest files in the `inputs/` directory (line 940-947)
3. When processing outputs, parent paths are checked against `path2inst` (line 989)

The assertion fails when **all items in `todo` are outputs with unresolved parent dependencies**.

### Evidence from KofamScan Run (KAgM2QMS)

#### Input Manifest Contents
```
IiGuFCzd → /msm_home/.../QUjpSrwaG84b/small_orfs.faa
CDoOod4a → /msm_home/.../V6CFBbwS8wr1/profiles
jh33jWJX → /msm_home/.../V6CFBbwS8wr1/ko_list
GjtEozCH → /msm_home/.../zQR7Bwyk3u52/kofamscan.oci
```

#### Library Manifest for QUjpSrwaG84b (`_metadata/index.yml`)
```yaml
manifest:
  ko_list:
    type: annotation::kofamscan_ko_list
  profiles:
    type: annotation::kofamscan_profiles
# NOTE: small_orfs.faa is NOT registered!
```

#### Physical Files in QUjpSrwaG84b Directory
```
ko_list  profiles  small_orfs.faa
```

### The Mismatch

| Source | small_orfs.faa Present? |
|--------|-------------------------|
| Input manifest (IiGuFCzd) | ✓ Points to `QUjpSrwaG84b/small_orfs.faa` |
| Physical file | ✓ Exists in directory |
| Library manifest | ✗ **NOT registered** |
| `path2inst` dict | ✗ **NOT included** (because `lib.Iterate()` only yields manifest items) |

### Failure Sequence

1. Given inputs (profiles, ko_list, kofamscan.oci) are processed first → added to `to_del`
2. Next iteration: only `kofamscan_results` output remains in `todo`
3. Output checks parents:
   - Parent `IiGuFCzd` (small_orfs.faa) → path from `kv2path`
   - Check: `ppath in path2inst`? → **NO** (not in library manifest)
   - Result: `ok = False`
4. Output is skipped (`continue`)
5. **`to_del` is empty** → assertion fails

## Proposed Fix

### Option 1: Add input manifest paths to `path2inst` (Recommended)

After building `kv2path` from input manifests, also add those paths to `path2inst`:

```python
# Around line 947, after building kv2path from input manifests:
for in_manifest in (output_path.parent/"inputs").iterdir():
    k = in_manifest.name
    with open(in_manifest) as f:
        for l in f:
            p = Path(l[:-1])
            _hash = md5(str(p).encode()).hexdigest()
            _hash = int(_hash[:15], 16)
            kv2path[(k, _hash)] = p, {}
            # FIX: Also add to path2inst if not already present
            if p not in path2inst and k in k2inst:
                path2inst[p] = k2inst[k]
```

### Option 2: Resolve path mismatch when checking parents

When checking if `ppath in path2inst`, also check with symlink resolution:

```python
# Around line 989:
ppath, _ = kv2path[k]
# Try both the literal path and resolved path
if ppath not in path2inst:
    ppath_resolved = ppath.resolve() if ppath.exists() else ppath
    if ppath_resolved not in path2inst:
        ok = False
        break
    ppath = ppath_resolved
```

### Option 3: Fix the staging to ensure library manifests match input manifests

The deeper issue is that files are staged (copied) but not properly registered in the library manifest. This would require changes to `StageWorkflow()` to ensure consistency.

## Additional Notes

- **DeepEC and ProteinBERT tests pass** because they don't use external binds - all inputs are staged consistently
- **KofamScan uses external binds** for profiles/ko_list directories, which may cause staging inconsistencies
- The output file (`kofamscan_results`) is valid and contains correct KO annotations - only the result compilation fails

## Files Affected

- `metasmith/agents.py` - `RunWorkflow()` function, lines ~940-1016
- Potentially `StageWorkflow()` for the root cause fix
