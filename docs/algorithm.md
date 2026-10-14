# Geometry workflow and preservation notes

## Original workflow

1. Read head vertices and facets; identify the skin using boundary marker `1005`.
2. Add concentric spherical tumor interfaces, marked `1007` and `1008`.
3. Resolve the selected location to a skin vertex and estimate its normal.
4. Construct the ceramic/electrode stack in a local tangent frame.
5. Find nearby skin triangles and build their edge-neighbor relationships.
6. Trace the electrode footprint through those triangles, split intersected
   facets, create the hydrogel side walls, and remove replaced faces.
7. Export the resulting PLC and region seed points for external TetGen meshing.

The geometry engine retains global node/facet dictionaries and the original
intersection logic. It is a research artifact rather than a general-purpose
mesh-processing library. The reader and renderer are separate, import-safe tools.

## Polygonal approximation of circular electrodes

`build_electrode_stack()` samples the circumference in the local tangent frame.
`div_c` includes a repeated closing point, so `append_electrode_geometry()` stores
`div_c - 1` distinct perimeter vertices per layer. Adjacent vertices form straight
edges: `div_c = 5` produces a square, while `div_c = 33` produces a 32-sided
approximation of a circle. Increasing this resolution refines the electrode
outline and adds geometry for the contact-boundary tracing and stitching steps.
The radius `a`, layer thicknesses, and array layout are separate parameters.

## Changes made during the refresh

- Consolidated the maintained editor from the repository's March 2021 script.
- Removed obsolete imports, commented-out experiments, and duplicate snapshots.
- Replaced comments/docstrings with English descriptions and normalized formatting.
- Moved initialization behind an explicit entry point; importing does not load data
  or open a figure.
- Added explicit input/output paths, gzip input, CLI commands, and installation metadata.
- Display the skin surface when opening the editor and explicitly enter the GUI loop.
- Fixed facet-length checks that compared the entire dictionary length rather than
  the selected facet, and made array-spacing selection repeat after four clicks.
- Set circumference sampling to 10, as in the later local experiment, instead of 4.
- Close output files explicitly.
- Use descriptive function names and separate the click timer from module imports.

The later local script contained an unfinished retry loop that could run forever
without recalculating the intersection. That draft was not adopted. The original
version remains available in Git history.

## Known limits

- GUI picking depends on Matplotlib's formatted coordinate text and the current
  view. It is not a robust ray/mesh intersection picker.
- Local neighbor searches are quadratic and may be slow on dense patches.
- Surface stitching has location-dependent failure cases. There is no rollback
  transaction for partial placement failure; restart the editor after an error.
- Initialize the editor only once per process. Mutable global state has not been
  comprehensively redesigned for concurrent sessions.
- The historical head reader assumes contiguous one-based IDs and triangular,
  single-polygon facets. Use the supplied marked input, not an arbitrary PLC.
- Output region seeds come from the original implementation. The base head has
  no region records; complete tissue assignments must be checked for a volume solve.
- The preview reader supports inline 3-D points and hole-free facet polygons.
  It triangulates polygons as fans for display, which assumes convex polygons.
  It is not a general TetGen validator or a watertightness/intersection test.
- PNG previews are sampled surface views. They do not demonstrate volumetric mesh
  quality, electrical accuracy, or clinical suitability.

## Verified model sizes

| Model | Points | Facets | Region records |
| --- | ---: | ---: | ---: |
| Base head (`mni4`) | 470,967 | 942,486 | 0 |
| Historical electrode example | 474,307 | 950,636 | 56 |

The input has boundary markers `1001` through `1006`. Only `1005` is explicitly
identified as skin by this code; other anatomical names are not inferred here.
The historical export also contains `-1`, `1`, `2`, `1007`, and `1008`.
