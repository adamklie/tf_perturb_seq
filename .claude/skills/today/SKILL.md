---
name: today
description: Daily session startup — read TODAY.md, check GitHub project status, present current state and priorities
user-invocable: true
---

Start-of-day routine for the TFP3 project. Do the following:

1. Read `TODAY.md` from the project root. Present the most recent date entry.
2. Check the GitHub project board status: `module load git && gh project item-list 4 --owner adamklie --format json`
3. Check open issues: `module load git && gh issue list --repo adamklie/tf_perturb_seq --state open --limit 20`
4. Summarize:
   - What was accomplished last session
   - What's currently running or blocked
   - What the open issues and priorities are
5. Ask Adam what he wants to focus on today.
6. Once decided, update TODAY.md with today's date and planned tasks.
