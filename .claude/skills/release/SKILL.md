---
name: release
description: Release the next fixed_point128 version from its dev/X branch - date the changelog, review the README, gate on tests, warnings and the documentation build, squash merge into main as one commit, tag, publish the GitHub release, and open the next dev branch. Use when asked to release, cut or publish a version, or merge a dev branch into main.
---

# fixed_point128 release

Turns the `dev/X` branch (X is the four part version, e.g. `1.1.0.0`) into release X: one squashed
commit on `main`, an annotated tag `vX`, and a GitHub release whose notes are the CHANGELOG entry.
This is how v0.11.0.0, v0.12.0.0 and v1.0.0.0 were released; keep every release the same.

Invoking `/release` is the go-ahead for the outward facing steps (push, PR, merge, tag, release).
Stop and ask only when a gate fails in a way that needs a decision, when the version is ambiguous,
and for the number of the next version in the last step.

## 0. Preconditions

- On branch `dev/X`, working tree clean, `gh auth status` logged in.
- X is the version in `include/fp128_shared.h` (`FP128_VERSION_MAJOR/MINOR/PATCH/BUILD`); it must
  match the branch name. The CMake `project(... VERSION ...)` in `CMakeLists.txt` is its first three
  parts.

## 1. CHANGELOG.md

- `## [Unreleased]` becomes `## [X] - <today, YYYY-MM-DD>`.
- At the bottom, the `[Unreleased]: .../compare/v<previous>...HEAD` link becomes
  `[X]: https://github.com/ericgur/fixed_point128/releases/tag/vX`.
- Read the section as release notes: no references to "this changelog" or "the commit before",
  no duplicate bullets, measurements stated against a named version.

## 2. README.md

Review it for facts the release made stale, and check each one against the code before changing it:
the "current version is" line in the `fp128_shared.h` section, behaviour and accuracy claims, the
Visual Studio project list, the `tools/` description.

## 3. Release gate - tests, warnings, documentation

```powershell
.claude\skills\release\scripts\Invoke-Fp128ReleaseCheck.ps1
```

Must end with `RELEASE CHECK PASSED`. It checks, and fails on any of:

| step | what passes |
| --- | --- |
| build | `fixed_point128_tests` and `fixed_point128_env_tests` built **from clean** with the `msvc` preset (Release and Debug) and `clang-cl` (Release), with **zero** compiler warnings |
| tests | both programs pass in all three builds |
| docs | `DoxyGen/Doxyfile` builds with **zero** warnings under the Doxygen version of `.github/workflows/docs.yml` - downloaded into `out\` if the local one differs - which is what the CI runs with warnings as errors |

The clean builds take 10-15 minutes. Logs are in `out\release-check\`. Fix whatever fails and run it
again; never release past a failure. `-SkipTests` / `-SkipDocs` rerun just one half while fixing.

## 4. Commit, push, pull request

1. Commit the CHANGELOG and README changes ("Prepare the X release notes and bring the README up to
   date"), push `dev/X`.
2. `gh pr create --base main --head dev/X --title "Dev/X merge" --body-file <file>`. Body: a one line
   headline, a numbered list of what the release brings (breaking changes named), a line stating
   the gate passed, then the footer `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.
3. Wait for the PR's checks: the `Documentation` workflow must pass (`gh pr checks <N>`), and
   `gh pr view <N> --json mergeStateStatus` must be `CLEAN`.

## 5. Squash merge into main

```sh
gh pr merge <N> --squash --subject "Dev/X merge (#<N>)" --body-file <file>
```

The body is the PR's numbered list and gate line, without the footer, ending with the
`Co-Authored-By` trailer. `main` gains exactly one commit. Keep the `dev/X` branch. Then confirm:
`git diff origin/main origin/dev/X` is empty, and the `Documentation` run on `main` for the merge
commit completed with success (it deploys the site).

## 6. Tag

```sh
git tag -a vX <merge commit> -F <message>   # then: git push origin vX
```

Message: `vX`, a blank line, a short wrapped summary that names the breaking changes, ending with
"See CHANGELOG.md for the full entry."

## 7. GitHub release

```powershell
.claude\skills\release\scripts\New-Fp128ReleaseNotes.ps1 -Version X -Ref vX
gh release create vX --verify-tag --title vX --latest --notes-file out\release-notes-X.md
```

The notes are the CHANGELOG section with `###` raised to `##`, wrapped lines joined and a link to the
full CHANGELOG at the tag - the same form as every earlier release. No assets.

## 8. Open the next development branch

Ask for the next version (default: the next minor, `M.(m+1).0.0`). Then:

1. `git switch -c dev/<next> origin/main`
2. Commit "Bump version to <next>": the four `FP128_VERSION_*` macros in `include/fp128_shared.h`.
3. Commit "Set the CMake project version to <M.m.p>": `project(... VERSION ...)` in `CMakeLists.txt`.
4. Commit "Open the Unreleased changelog section for <next>": `## [Unreleased]` above the release
   just made, and `[Unreleased]: https://github.com/ericgur/fixed_point128/compare/vX...HEAD` above
   its link.
5. `git push -u origin dev/<next>`

## Report

Finish with: the merge commit, the tag, the release URL, the gate results, and the new branch.
