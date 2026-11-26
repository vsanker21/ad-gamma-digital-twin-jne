# Zenodo supersession instructions (required manual step)

GitHub `main` and release `v2.0.0` are already corrected.  
Zenodo record **10.5281/zenodo.22746385** still serves the **flawed v1.0.0 zip** and must be superseded.

## Option A (recommended): New version of the existing concept DOI

1. Open https://zenodo.org/doi/10.5281/zenodo.22746385 while logged in as the depositor.
2. Click **New version**.
3. **Remove** the old `ad-gamma-digital-twin-jne_v1.0.0.zip` file from the new draft.
4. Upload a zip of `corrected_verified_deposit_v2` (or the GitHub `v2.0.0` source zip).
5. Update metadata:
   - Title: `Corrected verified archive ... (v2.0.0)`
   - Version: `2.0.0`
   - Description: state clearly that **v1.0.0 must not be used** (provenance / inclusion-rule / claim-boundary corrections).
6. Publish the new version. Concept DOI `10.5281/zenodo.22746384` will resolve to the latest version.

## Option B: API upload (if you have a Zenodo personal access token)

```powershell
$env:ZENODO_TOKEN = "<token>"
python "G:\Alzheimers Dementia - Digital Twin\Journal of Neural Engineering\corrected_verified_deposit_v2\code\upload_zenodo_new_version.py"
```

## After publish

Update `CITATION.cff` / README with the new version DOI (version DOI changes; concept DOI stays).
