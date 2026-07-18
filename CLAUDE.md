# CLAUDE.md

Claude Code entry point for this repo. The canonical, tool-agnostic brief is
**AGENTS.md** — imported below. Read it: it carries the working norms (critical
reviewer / statistician; verify don't assert), the two-phase pipeline invariants,
build/test commands, and where things live.

@AGENTS.md

## Claude-specific

- **Shared knowledge lives in [`.agents/`](.agents/)** — start at `.agents/README.md`.
  It has the roadmap, project state, the pigauto relationship, plans, and the
  simulation report. This directory is the durable, git-tracked record shared with
  the co-authors (Shinichi, Szymek).
- **Keep it current.** When you finish substantive work or learn something
  non-obvious, update `.agents/project-state.md` and `.agents/roadmap.md` so the next
  person or session inherits it — do not rely only on personal auto-memory.
- Personal auto-memory (`~/.claude/…/memory/`) is your working scratchpad; `.agents/`
  is the authoritative shared version.
