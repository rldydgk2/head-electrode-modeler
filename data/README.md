# Sample data

These files were supplied from the author's local research archive during the
2026 refresh. Compression is lossless; hashes below refer to the uncompressed
bytes. No vertex, facet, or region record was altered.

- `mni4.poly.gz`: the marked multi-surface input used by the historical editor.
- `electrode-example.poly.gz`: a saved March 2021 export containing head, tumor,
  and electrode geometry. It is included to make the README figures reproducible.

The older extraction script references `MNI152.mesh`, but the full derivation and
upstream license of `mni4.poly` are not recorded in the available source. No new
license is asserted for these anatomical meshes. Original TetGen binaries and
large volume-mesh outputs are not redistributed here.

| Distributed file | Original filename | Uncompressed SHA-256 |
| --- | --- | --- |
| `mni4.poly.gz` | `mni4.poly` | `3ce98d0d9037330a70f5115d017be855f4bd3d589e69e4a6e9b58e39cb74ce46` |
| `electrode-example.poly.gz` | `export3111027.poly` | `1aca214155230647146ac7a1fde2dd20e534e81837fbefb250ec6856207729b6` |
