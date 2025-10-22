# Report Generation - Folder Structure

This directory contains all documents related to the Advanced Mechanism Project Group 2 report generation.

## Folder Organization

### 📄 **01_Final_Documents**
- **ProjectReport1_Memo.pdf** - Final 6-page MEMO document (current version)
- **ProjectReport1_Updated.pdf** - Updated full report document

### 📝 **02_LaTeX_Source**
- **ProjectReport1_Memo.tex** - LaTeX source for the MEMO document
- **ProjectReport1_Updated.tex** - LaTeX source for the updated full report
- **references.bib** - Bibliography file with all citations

### 🖼️ **03_Images**
- **6 Bar options.png** - Six-bar mechanism topology options
- **Coupler Points.png** - Coupler point analysis for mechanism design
- **Example Watt 1.png** - Example Watt I mechanism configuration
- **Stephenson and Watt.png** - Stephenson and Watt chain configurations
- **Proposal Figure 1.jpg** - Proposal diagram 1
- **Proposal Figure 2.jpg** - Proposal diagram 2

### ⚙️ **04_Compilation_Files**
- All LaTeX auxiliary files (.aux, .bbl, .blg, .fdb_latexmk, .fls, .log, .out, .toc)
- These files are generated during LaTeX compilation and can be safely deleted if needed

### 📊 **05_Presentations**
- **Oral Presentation 1 AMD.pptx** - PowerPoint presentation

### 📦 **06_Archive**
- **ProjectReport1.pdf** - Original report document (archived)
- **ProjectReport1.tex** - Original LaTeX source (archived)

## Usage Notes

- **Current Working Document**: `01_Final_Documents/ProjectReport1_Memo.pdf`
- **Source Files**: All LaTeX source files are in `02_LaTeX_Source/`
- **Images**: All figures and diagrams are in `03_Images/`
- **Compilation**: Run LaTeX compilation from the `02_LaTeX_Source/` directory

## Compilation Instructions

To compile the MEMO document:
```bash
cd "02_LaTeX_Source"
latexmk -pdf -interaction=nonstopmode -halt-on-error "ProjectReport1_Memo.tex"
mv ProjectReport1_Memo.pdf "../01_Final_Documents/"
mv *.aux *.bbl *.blg *.fdb_latexmk *.fls *.log *.out "../04_Compilation_Files/"
```

To compile the updated report:
```bash
cd "02_LaTeX_Source"
latexmk -pdf -interaction=nonstopmode -halt-on-error "ProjectReport1_Updated.tex"
mv ProjectReport1_Updated.pdf "../01_Final_Documents/"
mv *.aux *.bbl *.blg *.fdb_latexmk *.fls *.log *.out *.toc "../04_Compilation_Files/"
```

**Note**: After compilation, always move the PDF to `01_Final_Documents/` and auxiliary files to `04_Compilation_Files/` to maintain the organized structure.
