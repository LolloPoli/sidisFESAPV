# sidisFESAPV
Sidis Italia shared code


# Basic workflow:


```
git clone https://github.com/nicolovalle/sidisFESAPV
```

Only once:
```
git config user.name "Your name"
git config user.email "your@email"
```

Edit code, then:
```
git add file.cpp
git commit -m "Your comment"
git push
```

## Content

### Macros
1. **epic_studies.cpp** runs over generated files (Pythia8), loops over all events, and fills a dedicated tree for the electron (MC and reco) and for the hadrons (MC and reco), storing all the relevant information (e.g. px, py, pz, PDG, Q2, xB, ... ). It takes as input a list file (.txt): the first line is the name of the output file, and each subsequent line is the path to an input file. To run the code from the terminal: 'root -l epic_studies.cpp\(\"lista_epic_studies.txt\"\)'.

2. **pion.plot2.cpp** which will become **plotting_relevant_plots.cpp** ...

### Workflow

To run it locally:
1. Move inside the environment ./eic-shell
2. Create the list file (.txt)
3. Run the epic_studies.cpp macro in order to store the information in the trees
4. Run the plotting_relevant_plots.cpp macro in order to ...

