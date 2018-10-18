# `firstpipe`

<!--
NOTE TO AUTHORS: replace all the **TODO** marks with your own content.
-->

**TODO**: insert high-level description of the pipeline.

The pipeline consists of the following stages:

  * [`example_run`](./example_run.sh). **TODO**: describe `example_run` stage.

  * [`tests`](./tests.sh). **TODO**: describe `tests` stage.

  * [`MS`](./MS.sh). **TODO**: describe `MS` stage.

# Obtaining the pipeline

To add this pipeline to your project using the
[`popper` CLI tool](https://github.com/systemslab/popper):

```bash
cd your-repo
popper add blubbundbla/physcraper/firstpipe
```

# Running the pipeline

To run the pipeline using the
[`popper` CLI tool](https://github.com/systemslab/popper):

```bash
cd physcraper
popper run firstpipe
```

The pipeline is executed on the following environment(s): `host`. In addition,
the following environment variables are expected:

  * `<ENV_VAR1>`. Description of variable.
  * `<ENV_VAR2>`. Another description.

> **TODO**: rename or remove ENV_VAR1 and ENV_VAR2 appropiately.

For example, the following is an execution with all expected
variables:

```bash
export <ENV_VAR1>=value-for-<ENV_VAR_1>
export <ENV_VAR2>=value-for-<ENV_VAR_2>

popper run firstpipe
```

> **TODO**: rename or remove `export` statements above appropriately.

# Dependencies

**TODO**: add list of dependencies, for example:

  * Python.
  * C++ compiler.
  * [Docker](https://docker.com) (for generating plots).
  * etc.
