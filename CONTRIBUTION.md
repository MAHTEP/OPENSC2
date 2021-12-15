# Contributig to OPENSC2


:fireworks::grinning:**Thanks for spending some time to contribute to the project**:grinning::fireworks:  

The following is a set of guidelines for contributing to [OPENSC2](https://github.com/MAHTEP/OPENSC2), which is hosted in the [MAHTEP Organization](https://github.com/MAHTEP) on GitHub. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request.

## How to Contribute?

To contribute to the project, please follow the following steps:

* Make a fork of the repositoty
* Clone your own forked repository
* Create a new brach
* Do your changes
* Test your changes
* Open a pull request following the instructions reported in the [template](https://github.com/MAHTEP/OPENSC2/blob/Contributions/.github/PULL_REQUEST_TEMPLATE/PULL_REQUEST_TEMPLATE.md)


## How to Open an Issue?

To open an issue:

* Clik the repository link
* Go to the [Issue](https://github.com/MAHTEP/OPENSC2/issues) section
* Search for you issue in the search toolbar
* If you do not find anthing related to your problem, please open a new issue selecting the most suitable template between _bug report_, _feature request_ and _custom template_. **Remember that blank issues are not allowed**.
* Fill the template
* Submit the issue

## <a name="commit"></a> Commit Message Format

We have very precise rules over how our Git commit messages must be formatted.
This format leads to **easier to read commit history**.

Each commit message consists of a **header**, a **body**, and a **footer**, as suggested by the [Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/)


```
<header>
<BLANK LINE>
<body>
<BLANK LINE>
<footer>
```

The `header` is mandatory and must conform to the [Commit Message Header](#commit-header) format.

The `body` is mandatory for all commits except for those of type "docs".
When the body is present it must be at least 20 characters long and must conform to the [Commit Message Body](#commit-body) format.

The `footer` is optional. The [Commit Message Footer](#commit-footer) format describes what the footer is used for and the structure it must have.


#### <a name="commit-header"></a>Commit Message Header

```
<type>(<scope>): <short summary>
  │       │             │
  │       │             └─⫸ Summary in present tense. Not capitalized. No period at the end.
  │       │
  │       └─⫸ Commit Scope: Prop_of_mat|Utility_functions|channel|conductors|coolant|environment|
  |                          fluid_components|jacket|mix_sc_stabilizer|opensc2_gui|simulation_starter|
  |                          simulations|solid_components|strands|super_conductor
  │
  └─⫸ Commit Type: build|docs|feat|fix|perf|refactor|test
```

The `<type>` and `<summary>` fields are mandatory, the `(<scope>)` field is optional.


##### Type

Must be one of the following:

* **build**: Changes that affect the build system or external dependencies
* **docs**: Documentation only changes
* **feat**: A new feature
* **fix**: A bug fix
* **perf**: A code change that improves performance
* **refactor**: A code change that neither fixes a bug nor adds a feature
* **test**: Adding missing tests or correcting existing tests


##### Scope
The scope should be the name of the affected module (as perceived by the person reading the changelog generated from commit messages).

The following is the list of supported scopes:

* `Prop_of_mat`
* `Prop_of_mat/aluminium`
* `Prop_of_mat/copper`
* `Prop_of_mat/epoxy`
* `Prop_of_mat/glass_epoxy`
* `Prop_of_mat/hastelloy_c276`
* `Prop_of_mat/jk2lb`
* `Prop_of_mat/incoloy_908`
* `Prop_of_mat/inconel_718`
* `Prop_of_mat/kapton`
* `Prop_of_mat/niobium_titanium`
* `Prop_of_mat/niobium3_tin`
* `Prop_of_mat/rare_earth_123`
* `Prop_of_mat/rare_earth_ba_cu_o`
* `Prop_of_mat/silver`
* `Prop_of_mat/sn60_pb40`
* `Prop_of_mat/stainless_steel`
* `Prop_of_mat/titanium`

* `Utility_functions`
* `Utility_functions\auxiliary_functions`
* `Utility_functions\gen_flow`
* `Utility_functions\initialization_functions`
* `Utility_functions\output`
* `Utility_functions\plots`
* `Utility_functions\solid_comp_initialization`
* `Utility_functions\transient_solution_functions`

* `channel`
* `conductors`
* `coolant`
* `environment`
* `fluid_components`
* `jacket`
* `mix_sc_stabilizer`
* `opensc2_gui`
* `simulation_starter`
* `simulations`
* `solid_components`
* `strands`
* `super_conductor`

##### Summary

Use the summary field to provide a succinct description of the change:

* use the imperative, present tense: "change" not "changed" nor "changes"
* don't capitalize the first letter
* no dot (.) at the end


#### <a name="commit-body"></a>Commit Message Body

Just as in the summary, use the imperative, present tense: "fix" not "fixed" nor "fixes".

Explain the motivation for the change in the commit message body. This commit message should explain _why_ you are making the change.
You can include a comparison of the previous behavior with the new behavior in order to illustrate the impact of the change.


#### <a name="commit-footer"></a>Commit Message Footer

The footer can contain information about breaking changes and deprecations and is also the place to reference GitHub issues, and other PRs that this commit closes or is related to.
For example:

```
BREAKING CHANGE: <breaking change summary>
<BLANK LINE>
<breaking change description + migration instructions>
<BLANK LINE>
<BLANK LINE>
Fixes #<issue number>
```
or

```
DEPRECATED: <what is deprecated>
<BLANK LINE>
<deprecation description + recommended update path>
<BLANK LINE>
<BLANK LINE>
Closes #<pr number>
```

Breaking Change section should start with the phrase "BREAKING CHANGE: " followed by a summary of the breaking change, a blank line, and a detailed description of the breaking change that also includes migration instructions.

Similarly, a Deprecation section should start with "DEPRECATED: " followed by a short description of what is deprecated, a blank line, and a detailed description of the deprecation that also mentions the recommended update path.


### Revert commits

If the commit reverts a previous commit, it should begin with `revert: `, followed by the header of the reverted commit.

The content of the commit message body should contain:

- information about the SHA of the commit being reverted in the following format: `This reverts commit <SHA>`,
- a clear description of the reason for reverting the commit message.



## How to Get Help?

If someting about contribution is not clear, please feel free to contact daniele.placido@polito.it. 

## Contributors License Agreement
Before sending pull request, you are encouraged to sign our [Contributors License Agreement (CLA)](https://github.com/MAHTEP/OPENSC2/blob/Contributions/CONTRIBUTOR_LICENSE_AGREEMENT.md). For any code changes to be accepted, the CLA must be signed. 
To do this:
- send an email to daniele.placido@polito.it to get the pdf version of the CLA
- read carefully and sign it
- send the signed document to laura.savoldi@polito.it (cc daniele.placido@polito.it)

If you have any question about CLA or the above procedure, write to daniele.placido@polito.it or laura.savoldi@polito.it.

<sub>The above procedure will be improved as soon as possible.<sub>

