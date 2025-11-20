
# 5. Troubleshooting

## 5.1. Restarting SaVor

If the workflow is killed or interrupted, append the `--rerun-incomplete` flag to the snakemake command. This allows regeneration of output files for rules that were in progress during the interruption.

```bash
snakemake -p --cores all --workflow-profile workflow-profiles/default --rerun-incomplete
```

If you had to manually interrupt the workflow or it was unexpectedly killed, the following command must first be executed before restarting SaVor:

```bash
snakemake -np --cores all --unlock
```