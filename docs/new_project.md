You might face the situation that you have to create a new `Sample_Project` within `spike`. Projects are essentially a means to organize your data and to control access to processing results via Snupy. Furthermore, a project can only contain **one** organism (human or mouse), might filter reads if being a xenograft, and produce coverage numbers for a set of gene panels and other quality control settings.

In general, a new project - let us call it `FooBar`, needs to be "registered" in the `config.yaml` file of `spike`. You will find the `projects:` section with existing projects like `Alps` or `Keimbahn`.
  - **project name**: Just add a new entry for `FooBar` at the end of the existing projects with the right indentation (i.e. one tab). This name must be precisely used in all sample sheets for samples belonging to this project, thus ensure it follows Illumina naming restrictions, e.g. no "-" or white spaces.
  - **species**: If we assume the project deals with human samples, you need to add the key-value pair `species: homo sapiens` "below", i.e. one more indent, the project name.
  - **snupy**: Since results of `spike` are uploaded to `snupy` and there is no common database, we here need to connect the `Sample_Project` (in spike) with the `project_id` (in snupy), which might be different for the two different instances running in `bonn` or at the `hhu`. Furthermore, one `contact` needs to be specified for every project in `snupy`, e.g. "Ute Fischer". You should now browse to the snupy instance of your choice and create an according new project `FooBar` there. Record the `project_id` (let us assume it is 42) and add a block to the config.yaml FooBar project like:
    ```
       snupy:
       hhu:
          contact: "Ute Fischer"
          project_id: 42
    ```
  - **genepanels**: The status update sheet can generate mean coverage values for given sets of genes. Those are defined as genepanels. Here, you can assign existing genepanels to the new project `FooBar` by adding the key `genepanels:` with the a list as value, e.g. 
    ```
       genepanels:
         - Sujal_Gosh
    ```
  - **min_coverage**: the above mentioned report file colours coverage values red if they fall below a certain threshold, e.g. 60. You can set this threshold here project specific, since we design libraries for human samples aiming for 60x coverage, mouse samples for 30x coverage.
  - **known_duos**: several projects include trio-computation, but sometimes we could not obtain samples from all three individuals. To mark those trios as incomplete, you can provide a list of `known_duos` to the project and list all `spike_entity_id` of incomplete trios.
  - **xenograft**: There are a few samples that were obtained by growing a xenograft in mice, e.g. human tissue was "implanted" into a mouse. In order to remove noise from the mouse background, reads get first mapped against the mouse genome. To activate this pre-processing step for project "FooBar" you can add the key-value pair `xenograft: "File.fa"` where `File.fa` is a hybrid multiple fasta file containing human and mouse genome data. ToDo: document how to create this file!
  - **capture_kit**: we are using two default capture kits for human and for mouse. This might change in the future and thus you can define a project specific capture_kit.

A complete example for the new project FooBar would look like this:

```
  FooBar:
    species: homo sapiens
    snupy:
      hhu:
        contact: "Ute Fischer"
        project_id: 42
    genepanels:
      - Sujal_Gosh
    min_coverage: 55
    known_duos:
      - family7
    xenograft: "GRCh_GRCm_filtering_hybrid.fa"
    capture_kit: SureSelect Human All Exon V6 r2
```
