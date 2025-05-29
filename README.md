# assignmant6

## Question 2 - Variant Code Robustness Comparison
### Goal
### Compare the robustness of real variant genetic codes to the Standard Genetic Code (SGC) by examining their local evolutionary landscapes. Tasks
1.	Hard-code two variant genetic codes:
    o	Vertebrate mitochondrial code
    o	CTG→Ser variant code in Candida
2.	Compute error costs:
    o	Apply the error-cost function described previously (Question 1) to these variant codes and the SGC.
3.	Simulate evolutionary neighborhood:
    o	Generate 1,000 random one-swap neighbors for each code (variant and SGC).
      1.	One-swap means to swap the designated amino acids of 2 codons
           o	Calculate error-costs for each neighbor.
4.	Visualization:
    o	Create boxplots to compare neighbor error-cost distributions for each genetic code.
5.	Analysis:
    o	Identify which genetic code resides on a "flatter" local optimum (i.e., more neighbor swaps have similar or improved robustness).
### Deliverables
•	Script/notebook with visualizations (boxplots and heatmap)
•	Brief report summarizing key insight
