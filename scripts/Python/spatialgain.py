# Packages
import numpy as np
import pandas as pd
import torch
import os
import re
import torch.nn.functional as F
from tqdm import tqdm

# Load data generated from R
# data_raw= pd.read_csv("../ml_io/data_preprocessed.csv").values
# mask_raw = pd.read_csv("../ml_io/mask.csv").values

print("SCRIPT STARTED")
print("Working dir:", os.getcwd())

input_dir = "../../output/SpatialGAIN"
output_dir = "../../output/results/SpatialGAIN"
os.makedirs(output_dir, exist_ok=True)

# Find all subdirectories that match the pattern (e.g., p10, p20, etc.)
subdirs = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d)) and re.match(r'p\d+', d)]

print("Detected folders:", subdirs)

# System Parameters
# 1. Mini batch size
mb_size = 128
# 2. Missing rate
# p_miss = 0.2
# 3. Hint rate
p_hint = 0.9
# 4. Loss Hyperparameters
alpha = 10
# 5. Train Rate
train_rate = 0.8

# Missing = mask_raw.astype(np.float32)
# Data = data_raw.astype(np.float32)

use_gpu = True  # set it to True to use GPU and False to use CPU

if use_gpu:
    torch.cuda.set_device(0)

for p_str in subdirs:
    print(f"\n{'='*30}\nProcessing missingness: {p_str}\n{'='*30}")
    
    # Construct paths based on the new structure
    current_path = os.path.join(input_dir, p_str)
    data_file_path = os.path.join(current_path, "data.csv")
    mask_file_path = os.path.join(current_path, "mask.csv")
    coords_file_path = os.path.join(current_path, "coords.csv")

    if not os.path.exists(coords_file_path):
        print(f"Skipping {p_str}: Missing coords.csv")
        continue

    if not os.path.exists(data_file_path) or not os.path.exists(mask_file_path):
        print(f"Skipping {p_str}: Missing data.csv or mask.csv")
        continue

    # 1. Load data and specific mask
    df_raw = pd.read_csv(data_file_path)
    data_raw = df_raw.values
    mask_raw = pd.read_csv(mask_file_path).values
    coords_raw = pd.read_csv(coords_file_path).values.astype(np.float32)
    # data_raw = pd.read_csv(os.path.join(input_dir, data_file_path)).values
    # mask_file = f"mask_{p_str}.csv"
    # mask_raw = pd.read_csv(os.path.join(input_dir, mask_file)).values

    Missing = mask_raw.astype(np.float32)
    Data = data_raw.astype(np.float32)

    # Spatial augmentation
    Coord_Dim = coords_raw.shape[1]  # should be 2 (x,y)

    Data_spatial = np.concatenate([Data, coords_raw], axis=1)

    # coordinates are always observed
    Coord_mask = np.ones((Missing.shape[0], Coord_Dim), dtype=np.float32)
    Missing_spatial = np.concatenate([Missing, Coord_mask], axis=1)

    Data = Data_spatial
    Missing = Missing_spatial

    No, Dim = Data.shape

    # Hidden state dimensions
    H_Dim1 = Dim
    H_Dim2 = Dim

    # Normalization (0 to 1) -> double normalisation but GAIN expects data to be in [0,1] range, so we need to normalise it first and then apply GAIN. After imputation, we can denormalise it back to the original scale.
    # Min_Val = np.zeros(Dim)
    # Max_Val = np.zeros(Dim)

    # for i in range(Dim):

    #     obs = Missing[:, i] == 1

    #     Min_Val[i] = np.min(Data[obs, i])
    #     Max_Val[i] = np.max(Data[obs, i])

    #     Data[:, i] = (Data[:, i] - Min_Val[i]) / (Max_Val[i] - Min_Val[i] + 1e-6)
    # Check if the data range is correctly normalised to [0,1]
    print("Data range:", Data.min(), Data.max())
    
    # Train Test Division     
    idx = np.random.permutation(No)

    Train_No = int(No * train_rate)
    Test_No = No - Train_No
        
    # Train / Test Features
    trainX = Data[idx[:Train_No],:]
    testX = Data[idx[Train_No:],:]

    # Train / Test Missing Indicators
    trainM = Missing[idx[:Train_No],:]
    testM = Missing[idx[Train_No:],:]

    # Necessary Functions

    # 1. Xavier Initialization Definition
    # def xavier_init(size):
    #     in_dim = size[0]
    #     xavier_stddev = 1. / tf.sqrt(in_dim / 2.)
    #     return tf.random_normal(shape = size, stddev = xavier_stddev)
    def xavier_init(size):
        in_dim = size[0]
        xavier_stddev = 1. / np.sqrt(in_dim / 2.)
        return np.random.normal(size = size, scale = xavier_stddev)
        

    # Hint Vector Generation
    def sample_M(m, n, p):
        A = np.random.uniform(0., 1., size = [m, n])
        B = A > p
        C = 1.*B
        return C

    # 1. Discriminator
    if use_gpu is True:
        D_W1 = torch.tensor(xavier_init([Dim*2, H_Dim1]),requires_grad=True, device="cuda")     # Data + Hint as inputs
        D_b1 = torch.tensor(np.zeros(shape = [H_Dim1]),requires_grad=True, device="cuda")

        D_W2 = torch.tensor(xavier_init([H_Dim1, H_Dim2]),requires_grad=True, device="cuda")
        D_b2 = torch.tensor(np.zeros(shape = [H_Dim2]),requires_grad=True, device="cuda")

        D_W3 = torch.tensor(xavier_init([H_Dim2, Dim]),requires_grad=True, device="cuda")
        D_b3 = torch.tensor(np.zeros(shape = [Dim]),requires_grad=True, device="cuda")       # Output is multi-variate
    else:
        D_W1 = torch.tensor(xavier_init([Dim*2, H_Dim1]),requires_grad=True)     # Data + Hint as inputs
        D_b1 = torch.tensor(np.zeros(shape = [H_Dim1]),requires_grad=True)

        D_W2 = torch.tensor(xavier_init([H_Dim1, H_Dim2]),requires_grad=True)
        D_b2 = torch.tensor(np.zeros(shape = [H_Dim2]),requires_grad=True)

        D_W3 = torch.tensor(xavier_init([H_Dim2, Dim]),requires_grad=True)
        D_b3 = torch.tensor(np.zeros(shape = [Dim]),requires_grad=True)       # Output is multi-variate

    theta_D = [D_W1, D_W2, D_W3, D_b1, D_b2, D_b3]

    # 2. Generator
    if use_gpu is True:
        G_W1 = torch.tensor(xavier_init([Dim*2, H_Dim1]),requires_grad=True, device="cuda")     # Data + Mask as inputs (Random Noises are in Missing Components)
        G_b1 = torch.tensor(np.zeros(shape = [H_Dim1]),requires_grad=True, device="cuda")

        G_W2 = torch.tensor(xavier_init([H_Dim1, H_Dim2]),requires_grad=True, device="cuda")
        G_b2 = torch.tensor(np.zeros(shape = [H_Dim2]),requires_grad=True, device="cuda")

        G_W3 = torch.tensor(xavier_init([H_Dim2, Dim]),requires_grad=True, device="cuda")
        G_b3 = torch.tensor(np.zeros(shape = [Dim]),requires_grad=True, device="cuda")
    else:
        G_W1 = torch.tensor(xavier_init([Dim*2, H_Dim1]),requires_grad=True)     # Data + Mask as inputs (Random Noises are in Missing Components)
        G_b1 = torch.tensor(np.zeros(shape = [H_Dim1]),requires_grad=True)

        G_W2 = torch.tensor(xavier_init([H_Dim1, H_Dim2]),requires_grad=True)
        G_b2 = torch.tensor(np.zeros(shape = [H_Dim2]),requires_grad=True)

        G_W3 = torch.tensor(xavier_init([H_Dim2, Dim]),requires_grad=True)
        G_b3 = torch.tensor(np.zeros(shape = [Dim]),requires_grad=True)

    theta_G = [G_W1, G_W2, G_W3, G_b1, G_b2, G_b3]


    # 1. Generator
    def generator(new_x,m):
        inputs = torch.cat(dim = 1, tensors = [new_x,m])  # Mask + Data Concatenate
        G_h1 = F.relu(torch.matmul(inputs, G_W1) + G_b1)
        G_h2 = F.relu(torch.matmul(G_h1, G_W2) + G_b2)   
        G_prob = torch.sigmoid(torch.matmul(G_h2, G_W3) + G_b3) # [0,1] normalized Output
        
        return G_prob


    # 2. Discriminator
    def discriminator(new_x, h):
        inputs = torch.cat(dim = 1, tensors = [new_x,h])  # Hint + Data Concatenate
        D_h1 = F.relu(torch.matmul(inputs, D_W1) + D_b1)  
        D_h2 = F.relu(torch.matmul(D_h1, D_W2) + D_b2)
        D_logit = torch.matmul(D_h2, D_W3) + D_b3
        D_prob = torch.sigmoid(D_logit)  # [0,1] Probability Output
        
        return D_prob

    # 3. Other functions
    # Random sample generator for Z
    def sample_Z(m, n):
        return np.random.uniform(0., 0.01, size = [m, n])        


    # Mini-batch generation
    def sample_idx(m, n):
        A = np.random.permutation(m)
        idx = A[:n]
        return idx


    def discriminator_loss(M, New_X, H):
        # Generator
        G_sample = generator(New_X,M)
        # Combine with original data
        Hat_New_X = New_X * M + G_sample * (1-M)

        # Discriminator
        D_prob = discriminator(Hat_New_X, H)

        # Loss
        D_loss = -torch.mean(M * torch.log(D_prob + 1e-8) + (1-M) * torch.log(1. - D_prob + 1e-8))
        return D_loss


    def generator_loss(X, M, New_X, H):
        # Structure
        # Generator
        G_sample = generator(New_X,M)

        # Combine with original data
        Hat_New_X = New_X * M + G_sample * (1-M)

        # Discriminator
        D_prob = discriminator(Hat_New_X, H)

        # Loss
        G_loss1 = -torch.mean((1-M) * torch.log(D_prob + 1e-8))
        MSE_train_loss = torch.mean((M * New_X - M * G_sample)**2) / torch.mean(M)

        G_loss = G_loss1 + alpha * MSE_train_loss 

        # MSE Performance metric
        MSE_test_loss = torch.mean(((1-M) * X - (1-M)*G_sample)**2) / torch.mean(1-M)
        return G_loss, MSE_train_loss, MSE_test_loss


    def test_loss(X, M, New_X):
        # Structure
        # Generator
        G_sample = generator(New_X,M)

        # MSE Performance metric
        MSE_test_loss = torch.mean(((1-M) * X - (1-M)*G_sample)**2) / torch.mean(1-M)
        return MSE_test_loss, G_sample

    optimizer_D = torch.optim.Adam(params=theta_D)
    optimizer_G = torch.optim.Adam(params=theta_G)

    # Start Iterations
    for it in tqdm(range(5000), desc=f"Training GAIN {p_str}"):    
        
        # Inputs
        mb_idx = sample_idx(Train_No, mb_size)
        X_mb = trainX[mb_idx,:]  
        
        Z_mb = sample_Z(mb_size, Dim) 
        M_mb = trainM[mb_idx,:]  
        H_mb1 = sample_M(mb_size, Dim, 1-p_hint)
        H_mb = M_mb * H_mb1
        
        New_X_mb = M_mb * X_mb + (1-M_mb) * Z_mb  # Missing Data Introduce
        
        if use_gpu is True:
            X_mb = torch.tensor(X_mb, device="cuda")
            M_mb = torch.tensor(M_mb, device="cuda")
            H_mb = torch.tensor(H_mb, device="cuda")
            New_X_mb = torch.tensor(New_X_mb, device="cuda")
        else:
            X_mb = torch.tensor(X_mb)
            M_mb = torch.tensor(M_mb)
            H_mb = torch.tensor(H_mb)
            New_X_mb = torch.tensor(New_X_mb)
        
        optimizer_D.zero_grad()
        D_loss_curr = discriminator_loss(M=M_mb, New_X=New_X_mb, H=H_mb)
        D_loss_curr.backward()
        optimizer_D.step()
        
        optimizer_G.zero_grad()
        G_loss_curr, MSE_train_loss_curr, MSE_test_loss_curr = generator_loss(X=X_mb, M=M_mb, New_X=New_X_mb, H=H_mb)
        G_loss_curr.backward()
        optimizer_G.step()    
            
    # Intermediate Losses
    if it % 100 == 0:
        print('Iter: {}'.format(it),end='\t')
        print('Train_loss: {:.4}'.format(np.sqrt(MSE_train_loss_curr.item())),end='\t')
        print('Test_loss: {:.4}'.format(np.sqrt(MSE_test_loss_curr.item())))

    # Testing (commented because we want to perform final imputation on the full dataset instead of just the test split)
    # Z_mb = sample_Z(Test_No, Dim) 
    # M_mb = testM
    # X_mb = testX
            
    # New_X_mb = M_mb * X_mb + (1-M_mb) * Z_mb  # Missing Data Introduce

    # if use_gpu is True:
    #     X_mb = torch.tensor(X_mb, device='cuda')
    #     M_mb = torch.tensor(M_mb, device='cuda')
    #     New_X_mb = torch.tensor(New_X_mb, device='cuda')
    # else:
    #     X_mb = torch.tensor(X_mb)
    #     M_mb = torch.tensor(M_mb)
    #     New_X_mb = torch.tensor(New_X_mb)
        
    # MSE_final, Sample = test_loss(X=X_mb, M=M_mb, New_X=New_X_mb)
            
    # print('Final Test RMSE: ' + str(np.sqrt(MSE_final.item())))

    # imputed_data = M_mb * X_mb + (1-M_mb) * Sample
    # print("Imputed test data:")
    # np.set_printoptions(formatter={'float': lambda x: "{0:0.8f}".format(x)})

    # if use_gpu is True:
    #     print(imputed_data.cpu().detach().numpy())
    # else:
    #     print(imputed_data.detach().numpy())

    # Final imputation on the full dataset

    # 1. Use the full dataset instead of just the test split
    Z_all = sample_Z(No, Dim)
    X_all = Data 
    M_all = Missing

    New_X_all = M_all * X_all + (1 - M_all) * Z_all

    # Move tensors to GPU if enabled
    if use_gpu:
        X_all_t = torch.tensor(X_all, device='cuda')
        M_all_t = torch.tensor(M_all, device='cuda')
        New_X_all_t = torch.tensor(New_X_all, device='cuda')
    else:
        X_all_t = torch.tensor(X_all)
        M_all_t = torch.tensor(M_all)
        New_X_all_t = torch.tensor(New_X_all)

    # 2. Run the Generator in "inference mode" (no gradient calculation)
    with torch.no_grad():
        G_sample_full = generator(New_X_all_t, M_all_t)

    # 3. Combine observed data with generated (imputed) data
    imputed_data_full = M_all_t * X_all_t + (1 - M_all_t) * G_sample_full
    imputed_np = imputed_data_full.detach().cpu().numpy()

    # remove coordinates from output
    imputed_np = imputed_np[:, :df_raw.shape[1]]

    # 4. Denormalise back to the R-preprocessed scale
    # for i in range(Dim):
    #     imputed_np[:, i] = imputed_np[:, i] * (Max_Val[i] - Min_Val[i] + 1e-6) + Min_Val[i]

    # 5. Restore Headers and Save
    # The columns from the original file are read to ensure the header is identical
    # original_headers = pd.read_csv("../ml_io/data_preprocessed.csv").columns
    # df_output = pd.DataFrame(imputed_np, columns=original_headers)
    df_res = pd.DataFrame(imputed_np, columns=df_raw.columns)
    output_filename = f"imputed_spatialgain_{p_str}.csv"
    df_res.to_csv(os.path.join(output_dir, output_filename), index=False)
    print(f"Saved: {output_filename}")
    # orig_headers = pd.read_csv(os.path.join(input_dir, data_file_path)).columns
    # df_res = pd.DataFrame(imputed_np, columns=orig_headers)

    # # df_output.to_csv("imputed_gain.csv", index=False)
    # # print(f"File saved! Dimensions: {df_output.shape}")
    # output_filename = f"imputed_gain_{p_str}.csv"
    # df_res.to_csv(os.path.join(output_dir, output_filename), index=False)
    # print(f"Saved: {output_filename} | Shape: {df_res.shape}")
