# Supplementary information

**OligoFormer: an accurate and robust prediction method for siRNA design**

Yilan Bai<sup>1,2</sup>, Haochen Zhong<sup>1,2</sup>, Taiwei Wang<sup>1,2,3</sup>, Zhi John Lu<sup>1,2, *</sup>

<sup>1</sup>MOE Key Laboratory of Bioinformatics, Center for Synthetic and Systems Biology, School of Life Sciences, Tsinghua University, Beijing, China.

<sup>2</sup>Institute for Precision Medicine, Tsinghua University, Beijing, China.

<sup>3</sup>Academy for Advanced Interdisciplinary Studies (AAIS), and Peking University–Tsinghua University–National Institute of Biological Sciences Joint Graduate Program (PTN), Peking University, Beijing, China

\*Corresponding author: Tel.: +86 10 62789217, E-mail: zhilu@tsinghua.edu.cn

## 1 Supplementary Methods

### 1.1 Data preprocessing

#### Redundancy Removal

Redundancy removal is a crucial step in constructing training and test sets to avoid the data leakage. We employed Needleman-Wunsch global alignment algorithm (Needleman and Wunsch 1970) to remove redundant RNA sequences between the training and test sets. The Needleman-Wunsch algorithm is a dynamic programming method widely used for sequence alignment, and allows for the optimal alignment of two siRNAs by maximizing the identity between them while considering possible gaps and mismatches. We calculated the identity percentage of each siRNA or mRNA sequence in the test set with every siRNA or mRNA sequence in the training set. We removed 1 siRNA from the Mixset that had more than 80% identity with any sequence in the Huesken set. Detailed results of redunancy are available on our Github library.

#### Efficacy normalization

The efficacy of Huesken, Mixset, and Takayuki datasets range from 0 to 1.341, -0.278 to 0.99, and 0 to 0.98, respectively. Min-max scaler was utilized to normalize the efficacy into a standardized inhibition efficiency scale ranging from 0 to 100%, which ensures that siRNA efficacy values were standardized as regression machine learning labels and comparable across datasets. We set the threshold at 70% inhibition efficiency for siRNAs in the dataset to classify siRNAs as positive or negative based on their inhibition efficiency.

#### Flanking region of target mRNAs

The selection of an appropriate flanking region surrounding the target mRNA sequence is important for efficient and specific mRNA embeddings. We systematically investigated the impact of flanking region length on siRNA efficacy prediction. We selected the complementary region of the 19 nt siRNA as the central point and extended it in both directions to obtain the corresponding mRNA sequences. This extension was conducted iteratively, ranging from 1 nt to 100 nt. Two hyperparameters “lm1” and “lm2” were employed to represent for the 5’ and 3’ flanking length to adapt to the different processing dimensions of the BiLSTM. We evaluated the siRNA efficacy prediction across different lengths of flanking region, and observed that flanking length of 19 nt outperformed the other symmetric flanking lengths, and in the analysis of asymmetric flanking length, the model stably performed well when both flanking regions were approximately 15-24 nt long (Supplementary Fig. 4). There are two possible reasons. Firstly, longer flanking regions provide more context of target mRNA but dilute the specific information of siRNA sequences, and vice versa. 19 nt flanking region may strike a balance between providing sufficient context of mRNA sequences and minimizing dilution of siRNA information. Secondly, the flanking length of 19 nt aligns with the processing length of siRNA, which facilitates the feature extraction of the BiLSTM. Therefore, 19 nt was chosen as the length of the flanking region of the mRNAs for the downstream analysis.

### 1.2 Thermodynamic parameters

The calculation of thermodynamic parameters is divided into two parts: frequency of single nucleotide or dinucleotide, and Gibbs free energy changes and enthalpy changes at specific position or over all siRNA. OligoFormer employs established algorithms to calculate thermodynamic parameters (Xia, SantaLucia et al. 1998, Ladunga 2007), including nearest-neighbor thermodynamic parameters, the Gibbs free energy (∆G), and enthalpy (∆H) of RNA duplex formation. The calculation equations for each thermodynamic parameter are described in detail below, and an example siRNA is shown for a better understanding of the calculation process.

#### Single nucleotide or dinucleotide

The single nucleotide N with the number k is whether the nucleotide at position k of siRNA is N, and the dinucleotide NM with the number k is whether the nucleotides at position k and position k+1 of siRNA are NM. Single nucleotide N with “all” is the frequency of nucleotide N in siRNA, and dinucleotide NM with “all” is the frequency of nucleotide NM in siRNA.

[Supplementary equation 1: single nucleotide and dinucleotide features N(k), NM(k), N(all), NM(all)](../assets/supp_figs/supplementary-equation-1.jpg)

$$\mathbf{N}(k)=\begin{cases}0, & siRNA[k]\neq\mathbf{N}\\ 1, & siRNA[k]=\mathbf{N}\end{cases},\quad \mathbf{N}\in[\mathbf{A},\mathbf{U},\mathbf{C},\mathbf{G}],\ k\in[1,19]$$

$$\mathbf{NM}(k)=\begin{cases}0, & siRNA[k:k+1]\neq\mathbf{NM}\\ 1, & siRNA[k:k+1]=\mathbf{NM}\end{cases},\quad \mathbf{N},\mathbf{M}\in[\mathbf{A},\mathbf{U},\mathbf{C},\mathbf{G}],\ k\in[1,19]$$

$$\mathbf{N}(all)=\frac{sum(\mathbf{N})}{19},\quad \mathbf{N}\in[\mathbf{A},\mathbf{U},\mathbf{C},\mathbf{G}]$$

$$\mathbf{NM}(all)=\frac{sum(\mathbf{NM})}{18},\quad \mathbf{N},\mathbf{M}\in[\mathbf{A},\mathbf{U},\mathbf{C},\mathbf{G}]$$

Take the first siRNA in Huesken dataset “CUAAUAUGUUAAUUGAUUU”,  for example,

[Supplementary equation 2: worked example of nucleotide and dinucleotide features](../assets/supp_figs/supplementary-equation-2.jpg)

$$\mathbf{U}(1)=0,\ \mathbf{G}(1)=0,\ \mathbf{UU}(1)=0,\ \mathbf{GG}(1)=0,\ \mathbf{U}(2)=1,\ \mathbf{C}(1)=1,\ \mathbf{CC}(1)=0,\ \mathbf{CG}(1)=0,\ \mathbf{A}(19)=0$$

$$\mathbf{U}(\mathrm{all})=\frac{10}{19}=0.526,\ \mathbf{G}(\mathrm{all})=\frac{2}{19}=0.105$$

$$\mathbf{GG}(\mathrm{all})=0,\ \mathbf{UA}(\mathrm{all})=\frac{3}{18}=0.167,\ \mathbf{CC}(\mathrm{all})=0,\ \mathbf{GC}(\mathrm{all})=0,\ \mathbf{UU}(\mathrm{all})=\frac{4}{18}=0.222$$

#### Gibbs free energy changes and enthalpy changes

Xia et al. (Xia, SantaLucia et al. 1998) calculated Gibbs free energy and enthalpy changes of all kinds of dinucleotide pairs for INN-HB Nearest-Neighbor Model(Supplementary Table 1.4). Gibbs free energy or enthalpy changes with the number k are the Gibbs free energy or enthalpy changes of the dinucleotide pair at the k and k+1 positions. While Gibbs free energy changes with “all” are a little bit more complicated to calculate using the following equation (Turner and Mathews 2010). The initiation term indicates siRNA-mRNA double strand structure formation, the per AU end is applied once per each AU pair at the end of siRNA, the symmetry term indicates self-complementary duplexes, and the last term is a sum of Gibbs free energy or enthalpy changes over all dinucleotide pairs.  The thermodynamic data needed to calculate these equations are shown in the Supplementary Table 1.4.

[Supplementary equation 3: ΔG°(k), ΔH°(k), ΔG°all, ΔH°all and ΔΔG°ends](../assets/supp_figs/supplementary-equation-3.jpg)

$$\Delta\mathbf{G}^{\circ}(k)=\Delta\mathbf{G}^{\circ}(siRNA[k:k+1]),\quad k\in[1,19]$$

$$\Delta\mathbf{H}^{\circ}(k)=\Delta\mathbf{H}^{\circ}(siRNA[k:k+1]),\quad k\in[1,19]$$

$$\Delta\mathbf{G}^{\circ}_{all}=\Delta\mathbf{G}^{\circ}_{initiation}+\Delta\mathbf{G}^{\circ}_{per\ AU\ end}+\Delta\mathbf{G}^{\circ}_{symmetry}+\sum_{k=1}^{18}\Delta\mathbf{G}^{\circ}(k)$$

$$\Delta\mathbf{H}^{\circ}_{all}=\Delta\mathbf{H}^{\circ}_{initiation}+\Delta\mathbf{H}^{\circ}_{per\ AU\ end}+\Delta\mathbf{H}^{\circ}_{symmetry}+\sum_{k=1}^{18}\Delta\mathbf{H}^{\circ}(k)$$

$$\Delta\Delta\mathbf{G}^{\circ}_{ends}=\Delta\mathbf{G}^{\circ}(siRNA[1:2])-\Delta\mathbf{G}^{\circ}(siRNA[18:19])+\Delta\mathbf{G}^{\circ}_{per\ AU\ end}$$

Take the first siRNA in Huesken dataset “CUAAUAUGUUAAUUGAUUU”,  for example,

[Supplementary equation 4: worked example of ΔG° and ΔH° values](../assets/supp_figs/supplementary-equation-4.jpg)

$$\Delta\mathbf{G}^{\circ}(1)=\Delta\mathbf{G}^{\circ}(siRNA[1:2])=\Delta\mathbf{G}^{\circ}(\mathbf{CU})=-2.08\ kcal/mol$$

$$\Delta\mathbf{H}^{\circ}(1)=\Delta\mathbf{H}^{\circ}(siRNA[1:2])=\Delta\mathbf{H}^{\circ}(\mathbf{CU})=-10.48\ kcal/mol$$

$$\Delta\mathbf{G}^{\circ}(2)=\Delta\mathbf{G}^{\circ}(siRNA[2:3])=\Delta\mathbf{G}^{\circ}(\mathbf{UA})=-1.33\ kcal/mol$$

$$\Delta\mathbf{G}^{\circ}(13)=\Delta\mathbf{G}^{\circ}(siRNA[13:14])=\Delta\mathbf{G}^{\circ}(\mathbf{UU})=-0.93\ kcal/mol$$

$$\Delta\mathbf{G}^{\circ}(18)=\Delta\mathbf{G}^{\circ}(siRNA[18:19])=\Delta\mathbf{G}^{\circ}(\mathbf{UU})=-0.93\ kcal/mol$$

$$\Delta\mathbf{H}^{\circ}_{all}=\Delta\mathbf{H}^{\circ}_{initiation}+\Delta\mathbf{H}^{\circ}_{per\ AU\ end}+\Delta\mathbf{H}^{\circ}_{symmetry}+\sum_{k=1}^{18}\Delta\mathbf{H}^{\circ}(k)$$

$$=3.61+3.72+0+(-10.48-7.69-6.82-9.38-7.69-9.38-10.44-11.40-6.82-7.69-6.82-9.38-6.82-10.44-12.44-9.38-6.82-6.82)$$

$$=-149.38\ kcal/mol$$

Therefore, the thermodynamic parameter vector of the example siRNA is :

[-1.6,-2.08,-10.48,0,0,-149.38,0.526,0,0.105,0,0,0,-1.33,0.167,1,1,0,-0.93,0,0,0,-0.93,0.222,0]

Please visit the Github library for more detailed codes.

### 1.3 RNA-FM embedding

RNA-FM (Zhang, Lang et al. 2024) (RNA Fundation Model) will transform RNA sequences into fixed-length feature vectors for the downstream feature fusion. The siRNA sequences and corresponding truncated mRNA sequences were extracted and formatted into fasta files, serving as input for RNA-FM module. We downloaded pre-trained models of RNA-FM from GoogleDrive (<https://drive.google.com/drive/folders/1VGye74GnNXbUMKx6QYYectZrY7G2pQ_J>), and applied embedding extraction function of RNA-FM to extract feature representations from the siRNA and mRNA sequences. The feature vectors were stored in npy (Numpy) file format, a binary file format used for storing numerical data in Python.  The shape of extracted RNA-FM embeddings is *l* × 640, where *l* is the length of siRNA or mRNA. Therefore, the shape of the siRNA embedding is 19 × 640 and the shape of the mRNA embedding is 57 × 640 for a flanking length of 19 nt. Next, Avg2D pooling was used to reduce the feature dimensions while preserving the most discriminative information. The kernel size of Avg2D pooling was set as *l* × 5, which is the most effective option, so that the 2D feature of *l* × 640 was transformed into the 1D feature of 128. This 1D vector of 128 is the RNA-FM embedding we need for the downstream feature fusion.

### 1.4 Oligo encoder embedding

The architecture of Oligo encoder consists of a 2D convolutional layer, a max pooling layer, an average pooling layer, a Bidirectional Long Short-Term Memory (BiLSTM) layer, two multi-head transformer encoder layers, and a flatten layer. Oligo encoder takes 19 nt siRNA sequences and 57 nt mRNA sequences **S** as input. We used one-hot encoding to represent the nucleotide sequences **S**, with each nucleotide (A, U, C, G, X) being assigned a unique binary vector of length 5. As a result, the one-hot embedding **E**<sub>O</sub> of one input sequence is a 1 × *l* × 5 matrix, where *l* is the length of siRNA or mRNA.

[Supplementary equation 5: one-hot embedding E_O](../assets/supp_figs/supplementary-equation-5.jpg)

$$\mathbf{E}_{O}=OneHot(\mathbf{S})=[[\mathbf{E}_{1},\mathbf{E}_{2},\mathbf{E}_{3},\cdots,\mathbf{E}_{l}]],\ \mathbf{E}_{i}\in\mathbb{R}^{5},\ i=1,2,3,\cdots,l$$

Next, we used a 2D convolutional layer to further extract features **E**<sub>**C**</sub>. The Conv2D layer was configured with an input channel of 1, an output channel of 64, a kernel size of 1 × 5 for siRNA and 5 × 5 for mRNA, a stride of 1, and no padding as described below:

[Supplementary equation 6: Conv2D feature extraction E_C](../assets/supp_figs/supplementary-equation-6.jpg)

$$\mathbf{E}_{\mathbf{C}}=Conv2D(\mathbf{E}_{O})=ReLU(\sum_{i=1}^{64}\mathbf{W}_{i}\times\mathbf{E}_{O})\in\mathbb{R}^{64\times l'\times 1}$$

$$(\mathbf{W}_{i}\times\mathbf{E}_{O})_{m,n}=\sum_{p=1}^{k_{1}}\sum_{q=1}^{k_{2}}\mathbf{W}_{i,p,q}\times{\mathbf{E}_{O}}_{m+p-1,n+q-1}$$

, where *l*′ represents 19 for siRNA and 53 for mRNA, (*k*<sub>1</sub>, *k*<sub>2</sub>) represents the kernel size, (**W**<sub>i</sub> × **E**<sub>O</sub>)<sub>m,n</sub> represents the value at position (*m*, *n*) in the output feature map, **W**<sub>i,p,q</sub> represents the weight at position (*p*, *q*) in the *i*-th kernel, and **E**<sub>O</sub><sub>m+p−1,n+q−1</sub> represents the value at position  (*m* + *p*−1, *n* + *q*−1) in the input feature map.

Then, a max pooling layer and an average pooling layer with a sliding window of fixed size 2 × 2 was applied to the transposed **E**<sub>C</sub>, which reduced the computational burden while retaining the most important information for subsequent layers in the network. The results of pooling layers **E**<sub>A</sub> ∈ ℝ<sup>l′×32</sup> and **E**<sub>C</sub> ∈ ℝ<sup>l′×32</sup> were concatenated with the original one-hot encoded feature **E**<sub>O</sub> ∈ ℝ<sup>l′×5</sup> to form a new feature matrix **E** ∈ ℝ<sup>l′×69</sup>.

[Supplementary equation 7: pooling and concatenation E_A, E_M, E](../assets/supp_figs/supplementary-equation-7.jpg)

$$\mathbf{E}_{A}=AvgPool_{2\times 2}(\mathbf{E}_{C}^{T})\in\mathbb{R}^{l'\times 32}$$

$$\mathbf{E}_{M}=MaxPool_{2\times 2}(\mathbf{E}_{C}^{T})\in\mathbb{R}^{l'\times 32}$$

$$\mathbf{E}=Concat(\mathbf{E}_{O},\mathbf{E}_{A},\mathbf{E}_{M})\in\mathbb{R}^{l'\times 69}$$

After the concatenation of **E**<sub>o</sub>, **E**<sub>A</sub> and **E**<sub>M</sub>, a 2-layer BiLSTM with the input dimension of 69 and the LSTM dimension of 32, was applied to **E**, described as follows:

[Supplementary equation 8: BiLSTM output E_B](../assets/supp_figs/supplementary-equation-8.jpg)

$$\mathbf{E}_{B}=dropout(ReLU(BiLSTM(\mathbf{E})))\in\mathbb{R}^{l'\times 64}$$

Multi-head self-attention (MHA) contained 8 heads, and each head was a mapping of query vector **Q**, key vector **K**, value vector **V**,  all of which were equal to **E**<sub>B</sub>. Suppose *N*<sub>h</sub> and *d*<sub>k</sub> denote the number of heads and the depth of keys, where *N*<sub>h</sub> = 8 and *d*<sub>k</sub> = 64. The output from each head *h* was computed as follows:

[Supplementary equation 9: attention head output A_h](../assets/supp_figs/supplementary-equation-9.jpg)

$$\mathbf{A}_{h}=dropout(softmax\left(\frac{\mathbf{E}_{B}\mathbf{W}_{q}(\mathbf{E}_{B}\mathbf{W}_{k})^{T}}{\sqrt{d_{k}}}\right)\mathbf{E}_{B}\mathbf{W}_{v},\quad h=1,2,\ldots,N_{h}$$

Then we concatenated the output from each head *h*. The  layer normalization (LN), feed-forward network (FFN), and residual connections were utilized to generate the output for each layer. We applied two linear layers with a ReLU activation layer in the FFN layer. Thus, the mechanism is summarized as:

[Supplementary equation 10: transformer encoder block O_0, z_n, O_n](../assets/supp_figs/supplementary-equation-10.jpg)

$$\mathbf{O}_{0}=[\mathbf{A}_{1},\mathbf{A}_{2},\ldots,\mathbf{A}_{N_{h}}])\in\mathbb{R}^{l'\times 64}$$

$$\mathbf{z}_{n}=MHA((LN(\mathbf{O}_{n-1})))+\mathbf{O}_{n-1}\in\mathbb{R}^{l'\times 64}$$

$$\mathbf{O}_{n}=FFN((LN(\mathbf{z}_{n})))+\mathbf{z}_{n}\in\mathbb{R}^{l'\times 64}$$

, where **O**<sub>n</sub> denotes the output from the *n*-th block and **z**<sub>n</sub> denotes the intermediate representation in the *n*-th block.

Finally **O**<sub>n</sub> was flattened into a vector **O** for the downstream feature fusion.

[Supplementary equation 11: flattened Oligo encoder output O](../assets/supp_figs/supplementary-equation-11.jpg)

$$\mathbf{O}=Flatten(\mathbf{O}_{n})\in\mathbb{R}^{l'*64}$$

### 1.5 Prediction calculation of OligoFormer

After the calculation module described above, we got thermodynamic embeddings ***T***, the RNA-FM embeddings **R**<sub>s</sub> and **R**<sub>m</sub> for the siRNA and mRNA sequences, and the embeddings **O**<sub>s</sub> and **O**<sub>m</sub> from the Oligo encoder for the siRNA and mRNA sequences, which was calculated as follows:

[Supplementary equation 12: module embeddings T, R_s/R_m, O_s/O_m](../assets/supp_figs/supplementary-equation-12.jpg)

$$\mathbf{T}=Thermodynamic(\mathbf{S}_{s},\mathbf{S}_{\mathbf{m}})\in\mathbb{R}^{24}$$

$$\mathbf{R}_{s},\mathbf{R}_{m}=RNAFM(\mathbf{S}_{s},\mathbf{S}_{\mathbf{m}})\in\mathbb{R}^{128}$$

$$\mathbf{O}_{s},\mathbf{O}_{m}=OligoEncoder(\mathbf{S}_{s},\mathbf{S}_{\mathbf{m}})\in\mathbb{R}^{l'*64}$$

, where *l*′=19 for siRNA and *l*′=53 for mRNA. The embeddings from these three modules were concatenated together to form a single feature vector **X**<sub>0</sub>.  The prediction calculation of OligoFormer involves a multilayer perceptron (MLP), which takes **X**<sub>0</sub> as input and output a scalar value *Ŷ* representing the predicted efficacy of the siRNA sequences.

[Supplementary equation 13: MLP prediction X_0, X_1, X_2, Ŷ](../assets/supp_figs/supplementary-equation-13.jpg)

$$\mathbf{X}_{0}=Concat(\mathbf{O}_{s},\mathbf{O}_{m},\mathbf{R}_{s},\mathbf{R}_{m},\mathbf{T})\in\mathbb{R}^{d_{0}}$$

$$\mathbf{X}_{1}=dropout(ReLU(\mathbf{X}_{0}\mathbf{W}_{1}+\mathbf{B}_{1}))\in\mathbb{R}^{d_{1}}$$

$$\mathbf{X}_{2}=dropout(ReLU(\mathbf{X}_{1}\mathbf{W}_{2}+\mathbf{B}_{2}))\in\mathbb{R}^{d_{2}}$$

$$\hat{Y}=softmax(\mathbf{X}_{2}\mathbf{W}_{3}+\mathbf{B}_{3})\in\mathbb{R}^{d_{3}}$$

, where **W**<sub>1</sub> ∈ ℝ<sup>d<sub>0</sub>×d<sub>1</sub></sup>,  **W**<sub>2</sub> ∈ ℝ<sup>d<sub>1</sub>×d<sub>2</sub></sup>,  **W**<sub>3</sub> ∈ ℝ<sup>d<sub>2</sub>×d<sub>3</sub></sup> are the weight matrices as parameters; *d*<sub>0</sub> = 4888 is the dimension of **X**<sub>0</sub> as MLP input; *d*<sub>1</sub> = 256,  *d*<sub>2</sub> = 64,  *d*<sub>3</sub> = 2 are numbers of neurons for three fully connected neural network layers; **B**<sub>1</sub> ∈ ℝ<sup>d<sub>1</sub></sup>, **B**<sub>2</sub> ∈ ℝ<sup>d<sub>2</sub></sup>, **B**<sub>3</sub> ∈ ℝ<sup>d<sub>3</sub></sup> are the terms of bias; ReLU is an activation function; softmax is a normalized exponential function; dropout is a dropout neural network layer with certain probability; *Ŷ* is the predicted siRNA efficacy.

### 1.6 Details of model training and test

#### Parameters of Model training

In the OligoFormer model, we utilized the PyTorch package in Python v3.7.16 to implement the network architecture. For model training and test, mean squared error (MSE) was chosen as the loss function, which quantified the average squared difference between predicted and actual efficacy values, penalizing larger errors more severely than smaller ones. We implemented Adam optimizer and weight decay function in the PyTorch package for optimized learning rate scheduling, which dynamically adjusted the learning rate during training and enhanced model convergence and performance. During model training, we optimized hyperparameters such as learning rate, weight decay, batch size, dropout probability, and epoch numbers inside the training set only. The learning rate was set to 0.0001, weight decay was set to 0.999 to ensure convergence, batch size was set to 16, the number of training epochs was set to 200 and early stopping was set to 30 epochs to prevent overfitting. Additionally, a fixed seed value of 42 was used to ensure reproducibility of results.

#### Evaluation metrics

Four evaluation metrics, area under the receiver operating characteristic curve (AUC), F1 score, area under the precision-recall curve (AUPRC), and Pearson correlation coefficient (PCC) were used to evaluate the prediction performance. By employing these training parameters and validation strategies, the OligoFormer model was trained effectively to accurately predict the efficacy of siRNA sequences. The following shows the calculation of these metrics for evaluation.

[Supplementary equation 14: MSE and Pearson correlation coefficient (PCC)](../assets/supp_figs/supplementary-equation-14.jpg)

$$MSE\ (Mean\ Square\ Error)=\frac{\sum_{i=1}^{n}(\hat{y}-y)^{2}}{n}$$

$$PCC=\frac{E(y\hat{y})-E(y)E(\hat{y})}{\sqrt{E(y^{2})-E^{2}(y)}\sqrt{E(\hat{y}^{2})-E^{2}(\hat{y})}}=\frac{\sum_{i=1}^{n}y\hat{y}-\frac{\sum_{i=1}^{n}y\sum_{i=1}^{n}\hat{y}}{n}}{\sqrt{(\sum_{i=1}^{n}y^{2}-\frac{(\sum_{i=1}^{n}y)^{2}}{n})(\sum_{i=1}^{n}\hat{y}^{2}-\frac{(\sum_{i=1}^{n}\hat{y})^{2}}{n})}}$$

[Supplementary equation 15: TP/TN/FP/FN definitions, precision and sensitivity (recall)](../assets/supp_figs/supplementary-equation-15.jpg)

$$TP\ (True\ Positive)=The\ number\ of\ correctly\ predicted\ positives$$

$$TN\ (True\ Negative)=The\ number\ of\ correctly\ predicted\ negatives$$

$$FP\ (False\ Positive)=The\ number\ of\ incorrectly\ predicted\ positives$$

$$FN\ (False\ Negative)=The\ number\ of\ incorrectly\ predicted\ negatives$$

$$Precision=\frac{TP}{TP+FP}$$

$$Sensitivity=Recall=\frac{TP}{TP+FN}$$

[Supplementary equation 16: specificity, FDR, accuracy and F1 score](../assets/supp_figs/supplementary-equation-16.jpg)

$$Specificity=\frac{TN}{TN+FP}$$

$$FDR=\frac{FP}{FP+TN}$$

$$Accuracy=\frac{TP+TN}{TP+FP+TN+FN}$$

$$F1\ score=2\times\frac{Precision\times Recall}{Precision+Recall}$$

#### Cross validation

During model test, we employed a 5-fold cross-validation strategy on the Huesken, Mixset, and Takayuki datasets for intra-dataset validation and we trained models on the Huesken dataset and evaluated the performance on the Mixset dataset for inter-dataset validation. We also performed cross cell line validation to ensure that our model could maintain high performance when applied to different cell lines, after excluding similar sequences. The Huesken dataset mainly included the H1299 cell line, with some other cell lines that did not coincide with cell lines of the other datasets and were called as the merged H1299 cell line (mH1299). Therefore, we focused on 5 different cell lines: mH1299, HaCaT, HEK293, HeLa, and T24. For each pair of cell lines, we trained the model on mH1299 dataset and tested it on another cell line dataset, ensuring the performance of OligoFormer was evaluated on unseen data from a different cell line. The results (Supplementary Table 1.5) demonstrate that OligoFormer maintains robust performance across different cell lines.

### 1.7 Comparison methods

For benchmarking, we compared six representative models for siRNA design, i.e., Monopoli-RF (Monopoli, Korkin et al. 2023), OligoWalk (Lu and Mathews 2008), siRNAPred (Kumar 2009), i-Score (Ichihara, Murakumo et al. 2007), s-Biopredsi (Huesken, Lange et al. 2005), and DSIR (Vert, Foveau et al. 2006). This section introduces how these comparison methods were implemented and the prediction results of test dataset are listed (Supplementary Table 3).

#### Monopoli-RF

Monopoli et al (Monopoli, Korkin et al. 2023). from Khvorova lab has recently published a machine learning approach to predict siRNA efficacy. This approach used asymmetric trichotomous grouping method with two independently selected thresholds: one threshold defined effective siRNAs and the other defined ineffective siRNAs, which split siRNA datsets into “effeactive”, “undefined”, and “ineffective” these three categories. Then it built a series of random forest (RF) models with k-fold cross-validation to form a final model to predict effective and ineffective siRNAs. The code for this model is available on Github (https://github.com/kmonopoli/trichotomous_data_partitioning_model_building). Since the author did not name the model, we called it Monopoli-RF for short. In order to maintain the consistency of the other comparison methods, the threshold for defining effective or ineffective siRNA  was set at 70% inhibition efficiency, thus there were no undefined siRNA partitions. We also trained the model on the original training data from the author's previous paper (Shmushkovich, Monopoli et al. 2018) and tested the performance (Supplementary Table 1.5).

#### OligoWalk

OligoWalk employs the nearest-neighbor thermodynamic model to calculate the binding free energy (∆G) of RNA duplex formation between an siRNA and the target mRNA. OligoWalk can predict the efficacy of siRNAs based on thermodynamic principles and molecular interactions. Because the OligoWalk software is embedded in the RNAstructure software, we downloaded RNAStructure (https://rna.urmc.rochester.edu/RNAstructure.html) and implemented OligoWalk model locally on the Linux platform. OligoWalk takes the fatsa file of the mRNA as input, and outputs efficacy predictions of all corresponding siRNAs. To reduce the calculation amount, we customed the start and end positions of mRNA by setting parameters --start and--end to the maximum range of siRNAs to be calculated. We set the length of predicted siRNA to 19, and obtained the efficacy predictions of siRNAs.

#### siRNAPred

siRNAPred involves support vector machines (SVMs) to predict efficacy of both 21mer and 19mer siRNAs. Nucleotide frequency features i.e. occurrence of mono to penta nucleotide string in a siRNA sequence were in develop hybrid-4 SVM method, binary patterns i.e. presence or absence of A, G, C, T/U nucleotides at each position in a siRNA sequence were used in binary pattern SVM method, and both methods were combined in hybrid-7 SVM method. We used the siRNAPred web server (https://webs.iiitd.edu.in/raghava/sirnapred/) to implement siRNAPred model. We chose efficacy prediction for 19mer and set the method as hybrid-7 method to get the efficacy predictions of siRNAs.

#### i-Score, s-Biopredsi, and DSIR

i-Score (inhibitory-Score) utilized a linear regression model on 2431 siRNAs to predict active siRNAs. It exclusively comprised of nucleotide preferences at each position and no other parameters. For s-Biopredsi, each siRNA sequence was input into the input layer neural network while the reporter data were utilized to update the weights between the network nodes with a learning rate of 0.1. This process was repeated ten times for each siRNA sequence and its corresponding target inhibition value. The final output was obtained by averaging the signals of the respective output nodes from all five networks. DSIR (Design of SIRna) is a simple linear model combining basic features of siRNA sequences for siRNA efficacy prediction. The i-Score Designer calculates nine different siRNA designing scores including i-Score, s-Biopredsi, and DSIR. And it provides an accessible web server (https://www.med.nagoya-u.ac.jp/neurogenetics/i_Score/i_score.html) where we can input mRNA and calculate efficacy predictions of 19 nt siRNA traversing the full length of mRNA online.

### 1.8 Ablation study

Ablation study was conducted to evaluate the impact of different features and modules on the performance of OligoFormer. It was divided into two parts: ablation study of different features, and ablation study of different modules (Supplementary Table 4).

#### Ablation study of different features

We examined the performance of different combinations of 4 types of features: siRNA sequences, mRNA sequences, RNA-FM embeddings, and thermodynamic parameters (TD). 15 feature sets were generated by combining these feature types: siRNA only, mRNA only, RNA-FM only, TD only, siRNA + mRNA, siRNA + RNA-FM, siRNA + TD, mRNA + RNA-FM, mRNA + TD, RNA-FM + TD, siRNA + mRNA + RNA-FM, siRNA + mRNA + TD, siRNA + RNA-FM + TD, mRNA + RNA-FM + TD, and combination of all 4 features. The results indicated that the full combination of all 4 features provided the best performance. The inclusion of RNA-FM features significantly enhanced the model's understanding of RNA sequences, while thermodynamic parameters contributed to the interpretation of molecular stability and interactions. A notable decline in performance was observed when features particularly RNA-FM and siRNA features were removed, which means these two features are crucial for capturing the complexities of RNA interactions. The thermodynamic parameters and mRNA features showed a less important impact when excluded, indicating that while they are beneficial, the model can still maintain a reasonable level of performance without them.

#### Ablation study of different modules

We focused on three basic modules of the OligoFormer architecture: Encoder layer, Conv2D, and BiLSTM. By creating models with different combinations of these modules, we evaluated the contribution of each to the overall performance, resulting in 7 module sets. Encoder layer Only, Conv2D only, BiLSTM only, Encoder layer + Conv2D, Encoder layer + BiLSTM, Conv2D + BiLSTM, and combination of all 3 modules. The full model achieved the best performance, demonstrating the complementary strengths of these modules. The Encoder layer as the core transformer encoder can capture the initial sequence representations, Conv2D can identify local patterns, and BiLSTM are crucial for capturing dependencies within the sequences. Removing any single component resulted in a drop in performance. Specifically, the absence of the Encoder layer led to the most significant decline, showing its critical role in initial sequence processing. The removal of Conv2D also significantly impacted performance. The impact of the exclusion of BiLSTM was slightly less severe compared to the other two modules.

### 1.9 Off-target

#### Functionality evaluation

To ensure the activity of the siRNA, some sequences with features known to be detrimental to RNAi functionality need to be removed (Birmingham, Anderson et al. 2007). (1) GC contents. High GC content may inhibit the dissociation of the duplex which is essential for RISC loading while very low GC content may cause lower target affinity and specificity. A range of 30–65% GC is considered optimal for identifying effective siRNAs and is generally utilized among potency-based siRNA design algorithms (Boudreau, Spengler et al. 2013). (2) five or more of any single bases in a row. Repeated bases have been shown to reduce functionality and selectivity and A or U stretches may additionally target regulatory motifs. (3) six or more consecutive G’s and/or C’s. Such regions have pronounced local stability which greatly inhibits duplex dissociation. (4) palindromic sequence. Internal complementary stretches more than three bases long may reduce the effective concentration of the mature duplex. Unstructured siRNAs result in higher efficacy. All of the above criteria are filtered by default.

#### Toxicity evaluation

Identification of putative toxic inducible factors within siRNAs provides a basis for the rational design of synthetic siRNAs that might avoid activating the innate immune response and improve the likelihood of passing toxicology tests. The toxicity of sequences are assessed and screened from two primary perspectives: (1) immunostimulatory motif. Toll-like receptors (TLRs) protect the host from pathogens by detecting infectious external agents including double-stranded RNA. Therefore, siRNAs may be recognized by the immune system as an unwanted invader and trigger an immune response (Jackson and Linsley 2010). Some known immunostimulatory motifs (GUCCUUCAA, UGUGU and CUGAAUU) will be filtered by default (Fakhr, Zare et al. 2016). (2) seed toxicity. Previous study has confirmed a cell type- and species- independent form of toxicity evoked by the 6mer seed sequence in siRNAs(Gao, Putzbach et al. 2018). All 4,096 6mer seed sequences and their corresponding cell viability have been downloaded (<https://www.6merdb.org/>). By default, siRNA 6mer seed regions with a cell viability of less than 50% are removed. Users can also customize the filter threshold ranging from 0 to 100%.

#### Incorporation of off-target search

PITA and TargetScan Context++ scores are utilized to evaluate the siRNA off-target effects. Both methods are capable of accepting 3'UTR sequences, siRNA sequences, and ORF sequences as inputs. UTR and ORF sequences are downloaded from TargetScan (<https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi>) . Note that the ORF sequences and the UTR sequences need to correspond with each other. Keep the default settings for other parameters.

If there is a possibility of siRNA and mRNA binding, the program will output a predicted off-target score (miRNA-like off-target score), otherwise it will not output a result, which means that this siRNA may not be off-target on this mRNA. For PITA, usually sequences with an off-target score below -10 have a greater chance of missing the target. For TargetsScan, there is no fixed threshold for determining whether the siRNA is off-target. This requires the users to choose the appropriate threshold based on practical needs and biological understanding. By default we set the threshold to 1. Users can also change the size of the two thresholds. Sequences that are predicted to miss the target by either method are filtered.

To test the off-target search module, off-target signatures of six different siRNAs targeting *APOB* gene in two human cell lines, HUH7 and PLC/PRF/5, were collected and defined as Burchard dataset (Burchard, Jackson et al. 2009). RNA was extracted at 6, 12, 24 and 48 hours upon transfection and gene expression was measured with microarray. We obtained the expression profile from the Gene Expression Omnibus database with accession number GSE14073 and performed the differential expression analysis through GEO2R interface. For each cell line and each siRNA, 6h and 12h samples were treated as technical replicates and log<sub>2</sub> fold changes (log<sub>2</sub>FC) of transcripts between mock and siRNA transfection were gained. Transcripts with log<sub>2</sub>FC below -0.4 are treated as off-target sequences and labeled as 1 while transcripts with log<sub>2</sub>FC between -0.4 and 0 were treated as off-target sequences and labeled as 0. We fed the sequences of these transcripts into the off-target prediction module and output the predictions with the default settings. Accuracy, sensitivity, specificity and F1-score of off-target prediction module were calculated.

### 1.10 Mismatch Handling Module

We provide an optional module to handle mismatches to enhance OligoFormer’s utility for designing siRNAs with differential activity on pathogenic alleles. This module is based on experimental data from this study (Becker, Ober-Reynolds et al. 2019), in which the dataset includes hsa-miR-21 and has-let-7a with specific mismatches and their corresponding association rate *k*<sub>on</sub>. And the ratio of *k*<sub>on</sub> of a mismatched sequence to that of the perfectly matched sequence serves as the label for our model. Each nucleotide sequence was one-hot encoded and converted into numerical matrices suitable for input of a neural network. A convolutional neural network (CNN) was employed to learn the impact of mismatches on binding efficiency, which took an original RNA and a mismatched RNA as input and output the *k*<sub>on</sub> ratio. More details of model training and test of mismatch module are on our Github repository. The final model was integrated into OligoFormer as an optional module. Users can select “-m” option and input a perfectly matched siRNA sequence along with a mismatched sequence of equal length. The module processes these inputs through the CNN and outputs the predicted *k*<sub>on</sub> ratio, indicating the binding efficiency of the mismatched sequence relative to the perfectly matched one.

### 1.11 User-friendly OligoFormer based on Docker

In order to create the docker image for OligoFormer, we made a customed dockerfile. We specified the base image as slim Python 3.8 and set the working directory to '/app' in the dockerfile, where OligoFormer and its dependencies will be installed. We then copied the OligoFormer source code and any required configuration files into the '/app' directory within the image. Next, we used the 'pip' package manager to install the necessary Python packages listed in OligoFormer's requirements.txt file. Additionally, we also installed some supplementary tools and dependencies required by OligoFormer to ensure its functionality within the docker container. Finally, we built the docker image using the 'docker build' command to complie the dockerfile instructions into a docker image encapsulating OligoFormer and its dependencies. We then uploaded the docker image of OligoFormer to the docker hub as a publicly accessible image.

With just two simple commands of pulling the OligoFormer docker image and running the OligoFormer container, users can quickly set up and utilize OligoFormer within a docker container environment. Firstly, fetch the latest OligoFormer docker image from the official docker hub repository by the command of “docker pull oligoformer:v1.0”. Secondly, create and run the docker container by the command of “docker run -it --name oligoformer-container -dt --restart unless-stopped oligoformer:v1.0 && docker exec -it oligoformer-container bash”. In this way, users can utilize OligoFormer conveniently on docker container environment in their workflows for oligonucleotide design.

## 2 Supplementary Figures

[Supplementary Fig. 1](../assets/supp_figs/supplementary-figure-1.jpg)

**Supplementary Fig. 1: Analysis of model loss and convergence speed.**

(a) Comparison of performance between original efficacy and normalized efficacy. (b) MSE loss of model for handling siRNA with some mismatches. (c) Loss change curves of models with different features. (d) AUC change curves of models with different features. The addition of mRNA features will slow down the convergence speed of the model, and the addition of RNA-FM and thermodynamic features will greatly accelerate the convergence of the model.

*Conversion note (not part of the source caption): in the published image, panel a shows the MSE loss for handling mismatches and panel b shows the original-versus-normalized efficacy comparison, i.e. the caption's (a)/(b) descriptions are swapped relative to the panels; the AUC change curves described as (d) are the unlabelled right-hand panel next to c.*

[Supplementary Fig. 2](../assets/supp_figs/supplementary-figure-2.jpg)

**Supplementary Fig. 2: AUC and PRC of different methods.**

AUC and PRC of Monopoli-RF, OligoWalk, siRNAPred, iScore, s-Biopredsi, DSIR, and OligoFormer, where (a) Huesken dataset, (b) Mixset, and (c) Takayuki dataset are intra-dataset validation, and (d) is inter-dataset validation.

[Supplementary Fig. 3](../assets/supp_figs/supplementary-figure-3.jpg)

**Supplementary Fig. 3: Performance of off-target prediction module.**

Accuracy (a), sensitivity (b), specificity (c), and F1-score (d) of off-target prediction module in Burchard dataset. For each of the six siRNAs from the Burchard dataset, we ran off-target prediction module for HUH7 and PLC/PRF/5 cell lines.

[Supplementary Fig. 4](../assets/supp_figs/supplementary-figure-4.jpg)

**Supplementary Fig. 4: Performance of symmetric and asymmetric flanking lengths.**

(a) AUC, PRC, and F1 score for symmetric flanking lengths range from 1 nt to 100 nt for inter-dataset validation. AUC (b), PRC (c), and F1 score (d) heatmap for asymmetric flanking lengths range from 1 nt to 100 nt for inter-dataset validation.

[Supplementary Fig. 5](../assets/supp_figs/supplementary-figure-5.jpg)

**Supplementary Fig. 5: Base preference maps of different datasets.**

Base preference maps for Huesken dataset (a), Mixset (b), and Takayuki dataset (c). siRNA is from 5 'to 3', mRNA is from 3 'to 5'. A larger absolute value of weight indicates greater importance for prediction, while a weight of zero indicates no importance of that base at that position.

[Supplementary Fig. 6](../assets/supp_figs/supplementary-figure-6.jpg)

**Supplementary Fig. 6: Generation and processing of siRNA efficacy.**

(a) Illustration of generation of siRNA efficacy, taking one target mRNA as an example. (b) Illustration of efficacy normalization during training and test, and inverse normalization during inference.

## References

Becker, W. R., B. Ober-Reynolds, K. Jouravleva, S. M. Jolly, P. D. Zamore and W. J. Greenleaf (2019). "High-Throughput Analysis Reveals Rules for Target RNA Binding and Cleavage by AGO2." Mol Cell **75**(4): 741-755 e711.

Birmingham, A., E. Anderson, K. Sullivan, A. Reynolds, Q. Boese, D. Leake, J. Karpilow and A. Khvorova (2007). "A protocol for designing siRNAs with high functionality and specificity." Nat Protoc **2**(9): 2068-2078.

Boudreau, R. L., R. M. Spengler, R. H. Hylock, B. J. Kusenda, H. A. Davis, D. A. Eichmann and B. L. Davidson (2013). "siSPOTR: a tool for designing highly specific and potent siRNAs for human and mouse." Nucleic Acids Res **41**(1): e9.

Burchard, J., A. L. Jackson, V. Malkov, R. H. Needham, Y. Tan, S. R. Bartz, H. Dai, A. B. Sachs and P. S. Linsley (2009). "MicroRNA-like off-target transcript regulation by siRNAs is species specific." RNA **15**(2): 308-315.

Fakhr, E., F. Zare and L. Teimoori-Toolabi (2016). "Precise and efficient siRNA design: a key point in competent gene silencing." Cancer Gene Ther **23**(4): 73-82.

Gao, Q. Q., W. E. Putzbach, A. E. Murmann, S. Chen, A. A. Sarshad, J. M. Peter, E. T. Bartom, M. Hafner and M. E. Peter (2018). "6mer seed toxicity in tumor suppressive microRNAs." Nat Commun **9**(1): 4504.

Huesken, D., J. Lange, C. Mickanin, J. Weiler, F. Asselbergs, J. Warner, B. Meloon, S. Engel, A. Rosenberg, D. Cohen, M. Labow, M. Reinhardt, F. Natt and J. Hall (2005). "Design of a genome-wide siRNA library using an artificial neural network." Nat Biotechnol **23**(8): 995-1001.

Ichihara, M., Y. Murakumo, A. Masuda, T. Matsuura, N. Asai, M. Jijiwa, M. Ishida, J. Shinmi, H. Yatsuya, S. Qiao, M. Takahashi and K. Ohno (2007). "Thermodynamic instability of siRNA duplex is a prerequisite for dependable prediction of siRNA activities." Nucleic Acids Research **35**(18).

Jackson, A. L. and P. S. Linsley (2010). "Recognizing and avoiding siRNA off-target effects for target identification and therapeutic application." Nat Rev Drug Discov **9**(1): 57-67.

Kumar, M., S. Lata, and G. P. S. Raghava (2009). "siRNApred: SVM based method for predicting efficacy value of siRNA." CSIR-IMTECH.

Ladunga, I. (2007). "More complete gene silencing by fewer siRNAs: transparent optimized design and biophysical signature." Nucleic Acids Res **35**(2): 433-440.

Lu, Z. J. and D. H. Mathews (2008). "OligoWalk: an online siRNA design tool utilizing hybridization thermodynamics." Nucleic Acids Res **36**(Web Server issue): W104-108.

Monopoli, K. R., D. Korkin and A. Khvorova (2023). "Asymmetric trichotomous partitioning overcomes dataset limitations in building machine learning models for predicting siRNA efficacy." Mol Ther Nucleic Acids **33**: 93-109.

Needleman, S. B. and C. D. Wunsch (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins." J Mol Biol **48**(3): 443-453.

Shmushkovich, T., K. R. Monopoli, D. Homsy, D. Leyfer, M. Betancur-Boissel, A. Khvorova and A. D. Wolfson (2018). "Functional features defining the efficacy of cholesterol-conjugated, self-deliverable, chemically modified siRNAs." Nucleic Acids Res **46**(20): 10905-10916.

Turner, D. H. and D. H. Mathews (2010). "NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure." Nucleic Acids Res **38**(Database issue): D280-282.

Vert, J. P., N. Foveau, C. Lajaunie and Y. Vandenbrouck (2006). "An accurate and interpretable model for siRNA efficacy prediction." BMC Bioinformatics **7**: 520.

Xia, T., J. SantaLucia, Jr., M. E. Burkard, R. Kierzek, S. J. Schroeder, X. Jiao, C. Cox and D. H. Turner (1998). "Thermodynamic parameters for an expanded nearest-neighbor model for formation of RNA duplexes with Watson-Crick base pairs." Biochemistry **37**(42): 14719-14735.

Zhang, Y., M. Lang, J. Jiang, Z. Gao, F. Xu, T. Litfin, K. Chen, J. Singh, X. Huang, G. Song, Y. Tian, J. Zhan, J. Chen and Y. Zhou (2024). "Multiple sequence alignment-based RNA language model and its application to structural inference." Nucleic Acids Res **52**(1): e3.

### Supplementary Table 1.1. siRNA datasets



[Supplementary Table 1.1. siRNA datasets](../assets/supp_table/supplementary-table-1-1.csv)

Source: Supplemntary Table 1.xlsx; worksheet: 1.1 siRNA datasets.
Merged cells retain their top-left value and blank covered cells; use merge ranges in metadata to interpret headers.
Merged ranges (top-left value retained): A1:C1.

| Supplementary Table 1.1. Three datasets of siRNAs |  |  |
| --- | --- | --- |
| Dataset | Number of siRNAs | Number of mRNAs |
| Huesken | 2431 | 34 |
| Takayuki | 702 | 1 |
| Mixset | 473 | 18 |

### Supplementary Table 1.2. Thermodynamic features



[Supplementary Table 1.2. Thermodynamic features](../assets/supp_table/supplementary-table-1-2.csv)

Source: Supplemntary Table 1.xlsx; worksheet: 1.2 Themodynamic features.
Merged cells retain their top-left value and blank covered cells; use merge ranges in metadata to interpret headers.
Merged ranges (top-left value retained): A1:C1.

| Supplementary Table 1.2. 24 thermodynamic features of siRNAs |  |  |
| --- | --- | --- |
| Individual feature |  Position | PCC* |
| ∆∆G°ends | 1 versus 19 | 0.46650000000000003 |
| ∆G° | 1 | 0.43480000000000002 |
| ∆H° | 1 | 0.40689999999999998 |
| U | 1 | 0.3634 |
| G | 1 | -0.3024 |
| ∆H° | all | 0.28129999999999999 |
| U | all | 0.24640000000000001 |
| UU | 1 | 0.24149999999999999 |
| G | all | -0.21510000000000001 |
| GG | 1 | -0.19270000000000001 |
| GC | 1 | -0.1988 |
| GG | all | -0.1925 |
| ∆G° | 2 | 0.18909999999999999 |
| UA | all | 0.1759 |
| U | 2 | 0.17610000000000001 |
| C | 1 | -0.17829999999999999 |
| CC | all | -0.17069999999999999 |
| ∆G° | 18 | -0.15939999999999999 |
| CC | 1 | -0.16139999999999999 |
| GC | all | -0.1613 |
| CG | 1 | -0.15559999999999999 |
| ∆G° | 13 | 0.153 |
| UU | all | 0.15079999999999999 |
| A | 19 | -0.12820000000000001 |
| *PCC (Pearson correlation coefficient) for efficacy and each thermodynamic features on the Hu dataset |  |  |

### Supplementary Table 1.3. Parameter settings



[Supplementary Table 1.3. Parameter settings](../assets/supp_table/supplementary-table-1-3.csv)

Source: Supplemntary Table 1.xlsx; worksheet: 1.3 Parameter settings.
Merged cells retain their top-left value and blank covered cells; use merge ranges in metadata to interpret headers.
Merged ranges (top-left value retained): A3:A7, A8:A15, A1:C1.

| Supplementary Table 1.3. Parameter settings of model and training |  |  |
| --- | --- | --- |
|  | Module | Parameter |
| Model parameters | 2D Convolution Layer  | Kernel size = 1x5/5x5, Stride = 1, Filters = 32, Pad = 0 |
|  | Max Pooling Layer |  Pool size = 2x2 |
|  | Average Pooling Layer  | Pool size = 2x2 |
|  | Bidirectional LSTM | Hidden dimension = 32 |
|  | Multi-head Transformer Encoder  | Layers: 2 layers, Hidden dimension = 64 |
| Training parameters | Optimizer | Adam |
|  | Loss function | MSE |
|  | Learning Rate | 1E-4 |
|  | Weight decay | 0.999 |
|  | Batch Size | 16 |
|  | Epochs | 200 |
|  | Early stop | 30 |
|  | Seed | 42 |

### Supplementary Table 1.4. Calculation reference



[Supplementary Table 1.4. Calculation reference](../assets/supp_table/supplementary-table-1-4.csv)

Source: Supplemntary Table 1.xlsx; worksheet: 1.4 Calculation reference.
Merged cells retain their top-left value and blank covered cells; use merge ranges in metadata to interpret headers.
Merged ranges (top-left value retained): A1:D1, B31:B32, C31:C32, D31:D32, B33:B34, C33:C34, D33:D34, B27:B28, C27:C28, D27:D28, B29:B30, C29:C30, D29:D30, B23:B24, C23:C24, D23:D24, B25:B26, C25:C26, D25:D26, B19:B20, C19:C20, D19:D20, B21:B22, C21:C22, D21:D22, B15:B16, C15:C16, D15:D16, B17:B18, C17:C18, D17:D18, B11:B12, C11:C12, D11:D12, B13:B14, C13:C14, D13:D14, B7:B8, C7:C8, D7:D8, B9:B10, C9:C10, D9:D10, B3:B4, C3:C4, D3:D4, B5:B6, C5:C6, D5:D6.

| Supplementary Table 1.4. RNA Thermodynamic Parameters for INN-HB Nearest-Neighbor Model* |  |  |  |
| --- | --- | --- | --- |
| parameters | ∆G°(kcal/mol) |  ∆H°(kcal/mol) | ∆S°(eu) |
| 5′AA3′ | -0.93 (0.03) | -6.82 (0.79) | -19.0 (2.5) |
| 3′UU5′ |  |  |  |
| 5′UU3′ | -0.93 (0.03) | -6.82 (0.79) | -19.0 (2.5) |
| 3′AA5′ |  |  |  |
| 5′AU3′ | -1.10 (0.08) | -9.38 (1.68) | -26.7 (5.2) |
| 3′UA5′ |  |  |  |
| 5′UA3′ | -1.33 (0.09) | -7.69 (2.02) | -20.5 (6.3) |
| 3′AU5′ |  |  |  |
| 5′CU3′ | -2.08 (0.06) | -10.48 (1.24) | -27.1 (3.8) |
| 3′GA5′ |  |  |  |
| 5′AG3′ | -2.08 (0.06) | -10.48 (1.24) | -27.1 (3.8) |
| 3′UC5′ |  |  |  |
| 5′CA3′ | -2.11 (0.07) | -10.44 (1.28) | -26.9 (3.9) |
| 3′GU5′ |  |  |  |
| 5′UG3′ | -2.11 (0.07) | -10.44 (1.28) | -26.9 (3.9) |
| 3′AC5′ |  |  |  |
| 5′GU3′ | -2.24 (0.06) | -11.40 (1.23) | -29.5 (3.9) |
| 3′CA5′ |  |  |  |
| 5′AC3′ | -2.24 (0.06) | -11.40 (1.23) | -29.5 (3.9) |
| 3′UG5′ |  |  |  |
| 5′GA3′ | -2.35 (0.06) | -12.44 (1.20) | -32.5 (3.7) |
| 3′CU5′ |  |  |  |
| 5′UC3′ | -2.35 (0.06) | -12.44 (1.20) | -32.5 (3.7) |
| 3′AG5′ |  |  |  |
| 5′CG3′ | -2.36 (0.09) | -10.64 (1.65) | -26.7 (5.0) |
| 3′GC5′ |  |  |  |
| 5′GG3′ | -3.26 (0.07) | -13.39 (1.24) | -32.7 (3.8) |
| 3′CC5′ |  |  |  |
| 5′CC3′ | -3.26 (0.07) | -13.39 (1.24) | -32.7 (3.8) |
| 3′GG5′ |  |  |  |
| 5′GC3′ | -3.42 (0.08) | -14.88 (1.58) | -36.9 (4.9) |
| 3′CG5′ |  |  |  |
| initiation | 4.09 (0.22) | 3.61 (4.12) | -1.5 (12.7) |
| per AU end |  0.45 (0.04) | 3.72 (0.83) | 10.5 (2.6) |
| symmetry correction(self-complementary) | 0.43 | 0 | -1.4 |
| symmetry correction(non-self-complementary) | 0 | 0 | 0 |
| *Xia, T., et al., Thermodynamic parameters for an expanded nearest-neighbor model for formation of RNA duplexes with Watson-Crick base pairs. Biochemistry, 1998. 37(42): p. 14719-35. |  |  |  |

### Supplementary Table 1.5. Cross validation



[Supplementary Table 1.5. Cross validation](../assets/supp_table/supplementary-table-1-5.csv)

Source: Supplemntary Table 1.xlsx; worksheet: 1.5 Cross validation.
CSV retains raw numeric values, not Excel display formatting (including percentages, dates and zero-padding); inspect source formats when interpreting units.
Merged cells retain their top-left value and blank covered cells; use merge ranges in metadata to interpret headers.
Merged ranges (top-left value retained): A3:D3, A11:D11, A1:D1.
Two blocks: cross cell line validation (rows 3-8) and inter-dataset validation (rows 11-20). The siRNAPred AUC and F1 score cells are stored in the source workbook as text beginning with an invisible U+FEFF character (retained verbatim; strip it before numeric parsing). Numeric cells keep raw stored precision (e.g. 0.74780000000000002 = 0.7478).

| Supplementary Table 1.5. Cross validation |  |  |  |  |
| --- | --- | --- | --- | --- |
|  |  |  |  |  |
| Cross cell line validation |  |  |  |  |
| Cell lines | AUC | PRC | F1 score |  |
| HaCaT | 0.8 | 0.74780000000000002 | 0.79169999999999996 |  |
| HEK293 | 0.82450000000000001 | 0.83299999999999996 | 0.82579999999999998 |  |
| HeLa | 0.77190000000000003 | 0.63480000000000003 | 0.55669999999999997 |  |
| T24 | 0.82620000000000005 | 0.43390000000000001 | 0.23810000000000001 |  |
|  |  |  |  |  |
|  |  |  |  |  |
| Inter-dataset validation |  |  |  |  |
| Methods | AUC | PRC | F1 score | PCC |
| Monopoli-RF | 0.74390000000000001 | 0.73880000000000001 | 0.70020000000000004 | 0.46279999999999999 |
| Monopoli-RF-ori | 0.60850000000000004 | 0.59640000000000004 | 0 | 0.4279 |
| OligoWalk | 0.70799999999999996 | 0.69979999999999998 | 0.62580000000000002 | 0.55920000000000003 |
| siRNAPred | ﻿0.5782 | 0.61799999999999999 | ﻿0.5245 | 0.20619999999999999 |
| iScore | 0.74619999999999997 | 0.75649999999999995 | 0.21970000000000001 | 0.36870000000000003 |
| s-Biopredsi | 0.73740000000000006 | 0.75570000000000004 | 0.59360000000000002 | 0.38490000000000002 |
| DSIR | 0.74829999999999997 | 0.76370000000000005 | 0.69450000000000001 | 0.503 |
| OligoFormer | 0.81630000000000003 | 0.81430000000000002 | 0.76900000000000002 | 0.58789999999999998 |

### Supplementary Table 2 (symmetric flanking lengths): AUC, AUPRC and F1 for three runs and average



[Supplementary Table 2 (symmetric flanking lengths): AUC, AUPRC and F1 for three runs and average](../assets/supp_table/supplementary-table-2-symmetric.csv)

Source: Supplemntary Table 2. Table of performance of different flanking lengths.xlsx; worksheet: symmetric.
Formula cells use saved cached values, which may be stale; formulas are preserved in metadata and are never evaluated.
Merged cells retain their top-left value and blank covered cells; use merge ranges in metadata to interpret headers.
Merged ranges (top-left value retained): B1:E1, F1:I1, J1:M1.
Row 1 merged group headers: AUC (B-E), AUPRC (F-I), F1 score (J-M); row 2: flanking length, runs 1-3, average. Averages are cached AVERAGE() formula values. Values keep raw stored precision.

|  | AUC |  |  |  | AUPRC |  |  |  | F1 score |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| flanking length | 1 | 2 | 3 | average | 1 | 2 | 3 | average | 1 | 2 | 3 | average |
| 1 | 0.81130000000000002 | 0.78380000000000005 | 0.81120000000000003 | 0.80209999999999992 | 0.81559999999999999 | 0.76839999999999997 | 0.80579999999999996 | 0.79660000000000009 | 0.75800000000000001 | 0.73540000000000005 | 0.75819999999999999 | 0.75053333333333327 |
| 2 | 0.7913 | 0.79469999999999996 | 0.79920000000000002 | 0.79506666666666659 | 0.79569999999999996 | 0.79620000000000002 | 0.78779999999999994 | 0.79323333333333323 | 0.74860000000000004 | 0.75319999999999998 | 0.74539999999999995 | 0.74906666666666666 |
| 3 | 0.8075 | 0.79800000000000004 | 0.8054 | 0.80363333333333342 | 0.81800000000000006 | 0.79820000000000002 | 0.79079999999999995 | 0.80233333333333334 | 0.75860000000000005 | 0.75319999999999998 | 0.75219999999999998 | 0.75466666666666671 |
| 4 | 0.79949999999999999 | 0.80620000000000003 | 0.79900000000000004 | 0.80156666666666665 | 0.81620000000000004 | 0.81410000000000005 | 0.79610000000000003 | 0.80880000000000007 | 0.71819999999999995 | 0.75439999999999996 | 0.75380000000000003 | 0.74213333333333331 |
| 5 | 0.80190000000000006 | 0.79730000000000001 | 0.79610000000000003 | 0.79843333333333344 | 0.79920000000000002 | 0.79010000000000002 | 0.79149999999999998 | 0.79360000000000008 | 0.76019999999999999 | 0.75749999999999995 | 0.74250000000000005 | 0.75340000000000007 |
| 6 | 0.80220000000000002 | 0.79879999999999995 | 0.80049999999999999 | 0.80049999999999999 | 0.81320000000000003 | 0.80379999999999996 | 0.79479999999999995 | 0.80393333333333328 | 0.76259999999999994 |  | 0.76590000000000003 | 0.76424999999999998 |
| 7 | 0.80020000000000002 | 0.7954 | 0.7944 | 0.79666666666666675 | 0.80510000000000004 | 0.79239999999999999 | 0.78410000000000002 | 0.79386666666666672 | 0.75039999999999996 | 0.73860000000000003 | 0.73609999999999998 | 0.74169999999999991 |
| 8 | 0.80779999999999996 | 0.80640000000000001 | 0.80640000000000001 | 0.80686666666666662 | 0.8075 | 0.80859999999999999 | 0.7903 | 0.80213333333333325 | 0.76900000000000002 | 0.76329999999999998 | 0.76380000000000003 | 0.76536666666666664 |
| 9 | 0.79020000000000001 | 0.79979999999999996 | 0.79420000000000002 | 0.79473333333333329 | 0.78180000000000005 | 0.80069999999999997 | 0.7833 | 0.78860000000000008 | 0.73680000000000001 | 0.75690000000000002 | 0.75549999999999995 | 0.74973333333333336 |
| 10 | 0.80079999999999996 | 0.80279999999999996 | 0.80620000000000003 | 0.80326666666666657 | 0.80059999999999998 | 0.80459999999999998 | 0.79879999999999995 | 0.80133333333333334 | 0.74350000000000005 | 0.74470000000000003 | 0.75439999999999996 | 0.74753333333333327 |
| 11 | 0.7964 | 0.8095 | 0.80659999999999998 | 0.8041666666666667 | 0.79349999999999998 | 0.81240000000000001 | 0.79530000000000001 | 0.80040000000000011 | 0.75 | 0.74790000000000001 | 0.75529999999999997 | 0.75106666666666666 |
| 12 | 0.80879999999999996 | 0.8054 | 0.80549999999999999 | 0.80656666666666654 | 0.81820000000000004 | 0.80700000000000005 | 0.79990000000000006 | 0.80836666666666668 | 0.75319999999999998 | 0.74639999999999995 | 0.75819999999999999 | 0.75260000000000005 |
| 13 | 0.79959999999999998 | 0.80290000000000006 | 0.79979999999999996 | 0.80076666666666663 | 0.79620000000000002 | 0.80110000000000003 | 0.77990000000000004 | 0.7924000000000001 | 0.75919999999999999 | 0.76729999999999998 | 0.75990000000000002 | 0.76213333333333333 |
| 14 | 0.81189999999999996 | 0.80790000000000006 | 0.80859999999999999 | 0.80946666666666667 | 0.81030000000000002 | 0.81620000000000004 | 0.80100000000000005 | 0.8091666666666667 | 0.76239999999999997 | 0.76259999999999994 | 0.75700000000000001 | 0.76066666666666671 |
| 15 | 0.79190000000000005 | 0.80920000000000003 | 0.79810000000000003 | 0.79973333333333352 | 0.78200000000000003 | 0.80659999999999998 | 0.78490000000000004 | 0.79116666666666668 | 0.74299999999999999 | 0.75270000000000004 | 0.76039999999999996 | 0.75203333333333333 |
| 16 | 0.79410000000000003 | 0.79600000000000004 | 0.79530000000000001 | 0.79513333333333336 | 0.78070000000000006 | 0.79490000000000005 | 0.7792 | 0.78493333333333337 | 0.74280000000000002 | 0.74209999999999998 | 0.74950000000000006 | 0.74480000000000002 |
| 17 | 0.81169999999999998 | 0.80820000000000003 | 0.80800000000000005 | 0.80930000000000002 | 0.81530000000000002 | 0.80569999999999997 | 0.79930000000000001 | 0.80676666666666674 | 0.75519999999999998 | 0.74870000000000003 | 0.76359999999999995 | 0.75583333333333336 |
| 18 | 0.81240000000000001 | 0.80549999999999999 | 0.80890000000000006 | 0.80893333333333339 | 0.80969999999999998 | 0.79330000000000001 | 0.80469999999999997 | 0.80256666666666676 | 0.75760000000000005 | 0.75449999999999995 | 0.75219999999999998 | 0.7547666666666667 |
| 19 | 0.80449999999999999 | 0.81130000000000002 | 0.82369999999999999 | 0.8131666666666667 | 0.81020000000000003 | 0.80610000000000004 | 0.81669999999999998 | 0.81099999999999994 | 0.76949999999999996 | 0.76329999999999998 | 0.7621 | 0.76496666666666668 |
| 20 | 0.8105 | 0.81220000000000003 | 0.80940000000000001 | 0.81070000000000009 | 0.80910000000000004 | 0.79520000000000002 | 0.79610000000000003 | 0.80013333333333347 | 0.76400000000000001 | 0.76290000000000002 | 0.75539999999999996 | 0.76076666666666659 |
| 21 | 0.80569999999999997 | 0.80159999999999998 | 0.80210000000000004 | 0.80313333333333325 | 0.79620000000000002 | 0.80130000000000001 | 0.79679999999999995 | 0.79810000000000014 | 0.7621 | 0.75 | 0.75580000000000003 | 0.75596666666666668 |
| 22 | 0.80479999999999996 | 0.80279999999999996 | 0.80820000000000003 | 0.80526666666666669 | 0.8 | 0.80490000000000006 | 0.76280000000000003 | 0.78923333333333334 | 0.75539999999999996 | 0.76119999999999999 | 0.72909999999999997 | 0.7485666666666666 |
| 23 | 0.8004 | 0.81040000000000001 | 0.81410000000000005 | 0.80830000000000002 | 0.8024 | 0.82100000000000006 | 0.81379999999999997 | 0.81240000000000012 | 0.75839999999999996 | 0.74519999999999997 | 0.76890000000000003 | 0.75749999999999995 |
| 24 | 0.80869999999999997 | 0.80320000000000003 | 0.78200000000000003 | 0.7979666666666666 | 0.81010000000000004 | 0.79369999999999996 | 0.74299999999999999 | 0.78226666666666667 | 0.75129999999999997 | 0.76290000000000002 | 0.75180000000000002 | 0.7553333333333333 |
| 25 | 0.80300000000000005 | 0.81030000000000002 | 0.81030000000000002 | 0.80786666666666684 | 0.79920000000000002 | 0.80930000000000002 | 0.80720000000000003 | 0.80523333333333336 | 0.75170000000000003 | 0.75519999999999998 | 0.75749999999999995 | 0.75479999999999992 |
| 26 | 0.81369999999999998 | 0.79749999999999999 | 0.79830000000000001 | 0.8031666666666667 | 0.81079999999999997 | 0.79139999999999999 | 0.77839999999999998 | 0.79353333333333331 | 0.76759999999999995 | 0.75939999999999996 | 0.7591 | 0.76203333333333323 |
| 27 | 0.77849999999999997 | 0.80049999999999999 | 0.80900000000000005 | 0.79599999999999993 | 0.77880000000000005 | 0.78510000000000002 | 0.79449999999999998 | 0.78613333333333335 | 0.73540000000000005 | 0.7571 | 0.74909999999999999 | 0.74719999999999998 |
| 28 | 0.81489999999999996 | 0.80249999999999999 | 0.80449999999999999 | 0.80730000000000002 | 0.81589999999999996 | 0.80069999999999997 | 0.78710000000000002 | 0.80123333333333335 | 0.75439999999999996 | 0.75039999999999996 | 0.75900000000000001 | 0.75459999999999994 |
| 29 | 0.7954 | 0.80079999999999996 | 0.81220000000000003 | 0.80280000000000007 | 0.79779999999999995 | 0.80110000000000003 | 0.80369999999999997 | 0.80086666666666673 | 0.75839999999999996 | 0.75319999999999998 | 0.74870000000000003 | 0.75343333333333329 |
| 30 | 0.80600000000000005 | 0.80490000000000006 | 0.80740000000000001 | 0.80609999999999993 | 0.78859999999999997 | 0.79349999999999998 | 0.79210000000000003 | 0.79139999999999999 | 0.75390000000000001 | 0.77010000000000001 | 0.75449999999999995 | 0.75950000000000006 |
| 31 | 0.8085 | 0.80459999999999998 | 0.80810000000000004 | 0.8070666666666666 | 0.80820000000000003 | 0.79779999999999995 | 0.79890000000000005 | 0.80163333333333331 | 0.76700000000000002 | 0.74609999999999999 | 0.74070000000000003 | 0.75126666666666664 |
| 32 | 0.80969999999999998 | 0.79910000000000003 | 0.80840000000000001 | 0.80573333333333341 | 0.80900000000000005 | 0.80030000000000001 | 0.79969999999999997 | 0.80300000000000005 | 0.76370000000000005 | 0.76619999999999999 | 0.75849999999999995 | 0.76280000000000003 |
| 33 | 0.80020000000000002 | 0.80459999999999998 | 0.81850000000000001 | 0.80776666666666674 | 0.79820000000000002 | 0.80169999999999997 | 0.81389999999999996 | 0.80459999999999987 | 0.7601 | 0.75690000000000002 | 0.76590000000000003 | 0.76096666666666657 |
| 34 | 0.79459999999999997 | 0.80200000000000005 | 0.80790000000000006 | 0.80149999999999999 | 0.79420000000000002 | 0.79949999999999999 | 0.79549999999999998 | 0.79640000000000011 | 0.7581 | 0.74909999999999999 | 0.74790000000000001 | 0.75170000000000003 |
| 35 | 0.81420000000000003 | 0.80810000000000004 | 0.79869999999999997 | 0.80700000000000005 | 0.80710000000000004 | 0.81100000000000005 | 0.7792 | 0.79910000000000003 | 0.75990000000000002 | 0.76039999999999996 | 0.75490000000000002 | 0.75839999999999996 |
| 36 | 0.81120000000000003 | 0.81500000000000006 | 0.81459999999999999 | 0.8136000000000001 | 0.80620000000000003 | 0.81769999999999998 | 0.80649999999999999 | 0.81013333333333326 | 0.73740000000000006 | 0.75670000000000004 | 0.76570000000000005 | 0.75326666666666675 |
| 37 | 0.79979999999999996 | 0.79269999999999996 | 0.79800000000000004 | 0.79683333333333328 | 0.7923 | 0.79049999999999998 | 0.7913 | 0.79136666666666666 | 0.75560000000000005 | 0.76870000000000005 | 0.74350000000000005 | 0.75593333333333346 |
| 38 | 0.80730000000000002 | 0.79900000000000004 | 0.81 | 0.80543333333333333 | 0.80979999999999996 | 0.80059999999999998 | 0.80359999999999998 | 0.80466666666666653 | 0.75790000000000002 | 0.75219999999999998 | 0.75 | 0.75336666666666663 |
| 39 | 0.80449999999999999 | 0.80300000000000005 | 0.79990000000000006 | 0.80246666666666666 | 0.80640000000000001 | 0.80020000000000002 | 0.7913 | 0.79930000000000001 | 0.75839999999999996 | 0.75360000000000005 | 0.7571 | 0.75636666666666663 |
| 40 | 0.79759999999999998 | 0.80730000000000002 | 0.80659999999999998 | 0.8038333333333334 | 0.78339999999999999 | 0.81369999999999998 | 0.79700000000000004 | 0.79803333333333326 | 0.75949999999999995 | 0.74729999999999996 | 0.75129999999999997 | 0.75269999999999992 |
| 50 | 0.80310000000000004 | 0.81310000000000004 | 0.80179999999999996 | 0.80600000000000005 | 0.80220000000000002 | 0.80820000000000003 | 0.78 | 0.79680000000000006 | 0.7571 | 0.75839999999999996 | 0.76459999999999995 | 0.76003333333333334 |
| 60 | 0.7994 | 0.80149999999999999 | 0.77710000000000001 | 0.79266666666666674 | 0.81320000000000003 | 0.7974 | 0.75849999999999995 | 0.78969999999999996 | 0.74639999999999995 | 0.74909999999999999 | 0.73850000000000005 | 0.7446666666666667 |
| 70 | 0.77890000000000004 | 0.80090000000000006 | 0.8004 | 0.79340000000000011 | 0.77470000000000006 | 0.79930000000000001 | 0.79300000000000004 | 0.78900000000000003 | 0.73870000000000002 | 0.74519999999999997 | 0.74429999999999996 | 0.74273333333333336 |
| 80 | 0.79220000000000002 | 0.79969999999999997 | 0.80530000000000002 | 0.79906666666666659 | 0.79010000000000002 | 0.79210000000000003 | 0.79869999999999997 | 0.7936333333333333 | 0.75270000000000004 | 0.74950000000000006 | 0.75660000000000005 | 0.75293333333333345 |
| 90 | 0.79900000000000004 | 0.78280000000000005 | 0.79890000000000005 | 0.79356666666666664 | 0.80210000000000004 | 0.77459999999999996 | 0.78159999999999996 | 0.78609999999999991 | 0.75680000000000003 | 0.76449999999999996 | 0.75180000000000002 | 0.75770000000000015 |
| 100 | 0.79259999999999997 | 0.79110000000000003 | 0.7611 | 0.78159999999999996 | 0.79280000000000006 | 0.77039999999999997 | 0.74850000000000005 | 0.77056666666666673 | 0.74260000000000004 | 0.75539999999999996 | 0.73299999999999998 | 0.74366666666666659 |

### Supplementary Table 2 (asymmetric flanking lengths): AUC grid, rows 5' flanking length x columns 3' flanking length



[Supplementary Table 2 (asymmetric flanking lengths): AUC grid, rows 5' flanking length x columns 3' flanking length](../assets/supp_table/supplementary-table-2-asymmetric-auc.csv)

Source: Supplemntary Table 2. Table of performance of different flanking lengths.xlsx; worksheet:  asymmetric_AUC.
Corner cell A1 reads "5' flanking3' flanking" in the source (diagonal header): row labels are 5' flanking lengths, column labels are 3' flanking lengths (1-40, then 50-100 by 10). Values keep raw stored precision.



### Supplementary Table 2 (asymmetric flanking lengths): AUPRC grid, rows 5' flanking length x columns 3' flanking length



[Supplementary Table 2 (asymmetric flanking lengths): AUPRC grid, rows 5' flanking length x columns 3' flanking length](../assets/supp_table/supplementary-table-2-asymmetric-prc.csv)

Source: Supplemntary Table 2. Table of performance of different flanking lengths.xlsx; worksheet:  asymmetric_PRC.
Corner cell A1 reads "5' flanking3' flanking" in the source (diagonal header): row labels are 5' flanking lengths, column labels are 3' flanking lengths (1-40, then 50-100 by 10). Values keep raw stored precision.



### Supplementary Table 2 (asymmetric flanking lengths): F1 grid, rows 5' flanking length x columns 3' flanking length



[Supplementary Table 2 (asymmetric flanking lengths): F1 grid, rows 5' flanking length x columns 3' flanking length](../assets/supp_table/supplementary-table-2-asymmetric-f1.csv)

Source: Supplemntary Table 2. Table of performance of different flanking lengths.xlsx; worksheet:   asymmetric_F1.
Corner cell A1 reads "5' flanking3' flanking" in the source (diagonal header): row labels are 5' flanking lengths, column labels are 3' flanking lengths (1-40, then 50-100 by 10). Values keep raw stored precision.



### Supplementary Table 3 (Huesken): per-siRNA measured efficacy and scores of OligoWalk, siRNAPred, iScore, s-Biopredsi, DSIR and OligoFormer



[Supplementary Table 3 (Huesken): per-siRNA measured efficacy and scores of OligoWalk, siRNAPred, iScore, s-Biopredsi, DSIR and OligoFormer](../assets/supp_table/supplementary-table-3-huesken.csv)

Source: Supplemntary Table 3. Table of siRNA prediction of benchmark methods.xlsx; worksheet: Huesken.
One row per 19-nt siRNA (antisense/guide sequence); 'efficacy' is the measured label; other columns are each method's predicted score.



### Supplementary Table 3 (Mixset): per-siRNA measured efficacy and scores of OligoWalk, siRNAPred, iScore, s-Biopredsi, DSIR and OligoFormer



[Supplementary Table 3 (Mixset): per-siRNA measured efficacy and scores of OligoWalk, siRNAPred, iScore, s-Biopredsi, DSIR and OligoFormer](../assets/supp_table/supplementary-table-3-mixset.csv)

Source: Supplemntary Table 3. Table of siRNA prediction of benchmark methods.xlsx; worksheet: Mixset.
CSV retains raw numeric values, not Excel display formatting (including percentages, dates and zero-padding); inspect source formats when interpreting units.
One row per 19-nt siRNA (antisense/guide sequence); 'efficacy' is the measured label; other columns are each method's predicted score.



### Supplementary Table 3 (Takayuki): per-siRNA measured efficacy and scores of OligoWalk, siRNAPred, iScore, s-Biopredsi, DSIR and OligoFormer



[Supplementary Table 3 (Takayuki): per-siRNA measured efficacy and scores of OligoWalk, siRNAPred, iScore, s-Biopredsi, DSIR and OligoFormer](../assets/supp_table/supplementary-table-3-takayuki.csv)

Source: Supplemntary Table 3. Table of siRNA prediction of benchmark methods.xlsx; worksheet: Takayuki.
Formula cells use saved cached values, which may be stale; formulas are preserved in metadata and are never evaluated.
One row per 19-nt siRNA (antisense/guide sequence); 'efficacy' is the measured label; other columns are each method's predicted score.



### Supplementary Table 3 (Inter-dataset): per-siRNA measured efficacy and scores of OligoWalk, siRNAPred, iScore, s-Biopredsi, DSIR and OligoFormer



[Supplementary Table 3 (Inter-dataset): per-siRNA measured efficacy and scores of OligoWalk, siRNAPred, iScore, s-Biopredsi, DSIR and OligoFormer](../assets/supp_table/supplementary-table-3-inter-dataset.csv)

Source: Supplemntary Table 3. Table of siRNA prediction of benchmark methods.xlsx; worksheet: Inter-dataset.
One row per 19-nt siRNA (antisense/guide sequence); 'efficacy' is the measured label; other columns are each method's predicted score.



### Supplementary Table 4. Ablation of OligoFormer components (check marks indicate included components per configuration column)



[Supplementary Table 4. Ablation of OligoFormer components (check marks indicate included components per configuration column)](../assets/supp_table/supplementary-table-4-ablation.csv)

Source: Supplemntary Table 4. Table of ablation.xlsx; worksheet: sheet1.
CSV retains raw numeric values, not Excel display formatting (including percentages, dates and zero-padding); inspect source formats when interpreting units.
Merged cells retain their top-left value and blank covered cells; use merge ranges in metadata to interpret headers.
Merged ranges (top-left value retained): B1:P1, B13:H13.

| component | Choice |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| siRNA | √ |  |  |  | √ | √ | √ |  |  |  | √ | √ | √ |  | √ |
| mRNA |  | √ |  |  | √ |  |  | √ | √ |  | √ | √ |  | √ | √ |
| RNA-FM |  |  | √ |  |  | √ |  | √ |  | √ | √ |  | √ | √ | √ |
| thermodynamics |  |  |  | √ |  |  | √ |  | √ | √ |  | √ | √ | √ | √ |
| AUC | 0.78900000000000003 | 0.75 | 0.70030000000000003 | 0.77569999999999995 | 0.79630000000000001 | 0.79395000000000004 | 0.79520000000000002 | 0.77010000000000001 | 0.78400000000000003 | 0.7853 | 0.80845 | 0.80095000000000005 | 0.80245 | 0.79500000000000004 | 0.81530000000000002 |
| PRC | 0.78029999999999999 | 0.75170000000000003 | 0.67900000000000005 | 0.76380000000000003 | 0.79720000000000002 | 0.79569999999999996 | 0.79190000000000005 | 0.73180000000000001 | 0.77270000000000005 | 0.78380000000000005 | 0.80920000000000003 | 0.79239999999999999 | 0.80110000000000003 | 0.7954 | 0.81430000000000002 |
| F1 score | 0.74104999999999999 | 0.72130000000000005 | 0.68689999999999996 | 0.73650000000000004 | 0.74544999999999995 | 0.75105 | 0.75065000000000004 | 0.74570000000000003 | 0.72670000000000001 | 0.73299999999999998 | 0.76554999999999995 | 0.75390000000000001 | 0.75160000000000005 | 0.73899999999999999 | 0.76890000000000003 |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| component | Choice |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| Conv2D | √ |  |  | √ | √ |  | √ |  |  |  |  |  |  |  |  |
| BiLSTM |  | √ |  | √ |  | √ | √ |  |  |  |  |  |  |  |  |
| OligoEncoder |  |  | √ |  | √ | √ | √ |  |  |  |  |  |  |  |  |
| AUC | 0.78610000000000002 | 0.77370000000000005 | 0.79710000000000003 | 0.78105000000000002 | 0.7893 | 0.80064999999999997 | 0.81530000000000002 |  |  |  |  |  |  |  |  |
| PRC | 0.78680000000000005 | 0.77310000000000001 | 0.79630000000000001 | 0.77110000000000001 | 0.78710000000000002 | 0.79479999999999995 | 0.81430000000000002 |  |  |  |  |  |  |  |  |
| F1 score | 0.73229999999999995 | 0.73860000000000003 | 0.75139999999999996 | 0.73145000000000004 | 0.75649999999999995 | 0.75444999999999995 | 0.76890000000000003 |  |  |  |  |  |  |  |  |

## Conversion notes

- Supplementary Materials (Methods, Supplementary Figs. 1-6, References) of the OligoFormer paper, converted from the publisher's Word file. Supplementary Tables 1-4 cited in the text were supplied as separate workbooks; they are included in this package as CSV files listed below under the supplementary tables.
- Displayed equations are preserved as image crops (supplementary-equation-1..16) each followed by a LaTeX transcription checked against the crop and the source .docx; source typos and unbalanced parentheses are preserved.
- Word margin line numbers were omitted as page furniture.
