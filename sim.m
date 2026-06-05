% =========================================================================
% 微波爐模糊控制器 - 完整比較版 (包含 COG, MOM, mMOM, CA 及 Excel 匯出)
% =========================================================================

% 1. 建立 Mamdani 模糊推論系統
fis = mamfis('Name', 'MicrowaveController');

% --- 定義輸入變數 ---
% 溫度 x: [-4, 4]
fis = addInput(fis, [-4 4], 'Name', 'Temp_x');
fis = addMF(fis, 'Temp_x', 'trimf', [-4 -4 0], 'Name', 'Low');
fis = addMF(fis, 'Temp_x', 'trimf', [-4 0 4], 'Name', 'Medium');
fis = addMF(fis, 'Temp_x', 'trimf', [0 4 4], 'Name', 'High');

% 重量 y: [0, 2]
fis = addInput(fis, [0 2], 'Name', 'Weight_y');
fis = addMF(fis, 'Weight_y', 'trimf', [0 0 1], 'Name', 'Light');
fis = addMF(fis, 'Weight_y', 'trimf', [0 1 2], 'Name', 'Medium');
fis = addMF(fis, 'Weight_y', 'trimf', [1 2 2], 'Name', 'Heavy');

% --- 定義輸出變數 ---
% 時間 z: [0, 10]
fis = addOutput(fis, [0 10], 'Name', 'Time_z');
fis = addMF(fis, 'Time_z', 'trimf', [0 0 5], 'Name', 'Short');
fis = addMF(fis, 'Time_z', 'trimf', [0 5 10], 'Name', 'Medium');
fis = addMF(fis, 'Time_z', 'trimf', [5 10 10], 'Name', 'Long');

% 功率 w: [600, 1200]
fis = addOutput(fis, [600 1200], 'Name', 'Power_w');
fis = addMF(fis, 'Power_w', 'trimf', [600 600 900], 'Name', 'Low');
fis = addMF(fis, 'Power_w', 'trimf', [600 900 1200], 'Name', 'Medium');
fis = addMF(fis, 'Power_w', 'trimf', [900 1200 1200], 'Name', 'High');

% --- 定義模糊規則 ---
% [Temp, Weight, Time, Power, Weight, AND(1)]
ruleList = [
    1 3 3 3 1 1;  % R1: Low & Heavy -> Long & High
    1 2 2 3 1 1;  % R2: Low & Medium -> Medium & High
    1 1 1 3 1 1;  % R3: Low & Light -> Short & High
    2 3 3 2 1 1;  % R4: Medium & Heavy -> Long & Medium
    2 2 2 2 1 1;  % R5: Medium & Medium -> Medium & Medium
    2 1 1 2 1 1;  % R6: Medium & Light -> Short & Medium
    3 3 3 1 1 1;  % R7: High & Heavy -> Long & Low
    3 2 2 1 1 1;  % R8: High & Medium -> Medium & Low
    3 1 1 1 1 1;  % R9: High & Light -> Short & Low
];
fis = addRule(fis, ruleList);

% =========================================================================
% 2. 建立測試網格
% =========================================================================
x_vec = -4:0.:4;
y_vec = 0:0.1:2;
[X, Y] = meshgrid(x_vec, y_vec);
inputData = [X(:), Y(:)]; % 將網格攤平為 Nx2 矩陣

% =========================================================================
% 3. 計算四種去模糊化方法的結果
% =========================================================================

% [方法 1] COG (重心法)
fis.DefuzzificationMethod = 'centroid';
out_cog = evalfis(fis, inputData);

% [方法 2] MOM (最大平均法)
fis.DefuzzificationMethod = 'mom';
out_mom = evalfis(fis, inputData);

% [方法 3] Modified MOM (修改型最大平均法)
fis.DefuzzificationMethod = 'som';
out_som = evalfis(fis, inputData);
fis.DefuzzificationMethod = 'lom';
out_lom = evalfis(fis, inputData);
out_mmom = (out_som + out_lom) / 2;

% [方法 4] CA (中心平均法) 
% 直接使用內建的 trimf 函數，算出 861 個點在各個模糊集合的歸屬度
T_Low = trimf(inputData(:,1), [-4 -4 0]);
T_Med = trimf(inputData(:,1), [-4 0 4]);
T_High = trimf(inputData(:,1), [0 4 4]);

W_Light = trimf(inputData(:,2), [0 0 1]);
W_Med = trimf(inputData(:,2), [0 1 2]);
W_Heavy = trimf(inputData(:,2), [1 2 2]);

% 根據 9 條規則計算激發強度 (w_i = min(mu_x, mu_y))
w1 = min(T_Low, W_Heavy);  % R1: Low & Heavy
w2 = min(T_Low, W_Med);    % R2: Low & Medium
w3 = min(T_Low, W_Light);  % R3: Low & Light
w4 = min(T_Med, W_Heavy);  % R4: Medium & Heavy
w5 = min(T_Med, W_Med);    % R5: Medium & Medium
w6 = min(T_Med, W_Light);  % R6: Medium & Light
w7 = min(T_High, W_Heavy); % R7: High & Heavy
w8 = min(T_High, W_Med);   % R8: High & Medium
w9 = min(T_High, W_Light); % R9: High & Light

% 將 9 條規則的強度組合為 861x9 矩陣
W = [w1, w2, w3, w4, w5, w6, w7, w8, w9];

% 輸出中心點 (改為 1x9 橫向量，方便後續純矩陣相乘)
C_Z = [10, 5, 0, 10, 5, 0, 10, 5, 0];
C_W = [1200, 1200, 1200, 900, 900, 900, 600, 600, 600];

% 避免分母為 0 (沒有規則被觸發時的防呆機制)
sum_w = sum(W, 2);
sum_w(sum_w == 0) = eps;

% 執行矩陣相乘完成加權平均公式
out_ca_z = (W * C_Z') ./ sum_w;
out_ca_w = (W * C_W') ./ sum_w;
out_ca = [out_ca_z, out_ca_w];

% =========================================================================
% 4. 準備資料匯出至 Excel
% =========================================================================
filename = 'Fuzzy_Simulation_Results.xlsx';
disp(['準備匯出資料至 ', filename, ' ...']);

precision_bits = 4;

% 將座標軸數值加上單位，轉換為字串 Cell
x_headers = arrayfun(@(x) sprintf('%g°C', x), x_vec, 'UniformOutput', false);
y_headers = arrayfun(@(y) sprintf('%gkg', y), y_vec', 'UniformOutput', false);

% 使用 size(X, 1) 固定列數，並用 [] 讓 MATLAB 自動推算對應的行數
% 這樣寫能完美避開版本不相容或元素數量判斷的問題
resultsMap = {
    'COG_Time',   reshape(out_cog(:, 1), size(X, 1), []);
    'COG_Power',  reshape(out_cog(:, 2), size(X, 1), []);
    'MOM_Time',   reshape(out_mom(:, 1), size(X, 1), []);
    'MOM_Power',  reshape(out_mom(:, 2), size(X, 1), []);
    'mMOM_Time',  reshape(out_mmom(:, 1), size(X, 1), []);
    'mMOM_Power', reshape(out_mmom(:, 2), size(X, 1), []);
    'CA_Time',    reshape(out_ca(:, 1), size(X, 1), []);
    'CA_Power',   reshape(out_ca(:, 2), size(X, 1), []);
};

% 寫入 Excel (writecell 支援字串與數字混編)
for i = 1:size(resultsMap, 1)
    sheetName = resultsMap{i, 1};

    % 核心數據四捨五入，並轉為 Cell 型態
    dataMatrix = round(resultsMap{i, 2}, precision_bits);
    dataMatrixSwapped = dataMatrix.';
    dataCell = num2cell(dataMatrixSwapped);
    
    % 組合最終輸出的 Cell 陣列
    % 第一列：[空白, 溫度標頭...]
    % 下方列：[重量標頭, 數據1, 數據2...]
    topRow = [{' '}, y_headers.'];
    dataRows = [x_headers.', dataCell];
    exportCell = [topRow; dataRows];
    
    % 匯出至指定的 Sheet
    writecell(exportCell, filename, 'Sheet', sheetName);
    fprintf('  - %s 寫入完成 (已附加單位)\n', sheetName);
end

disp('所有資料匯出完畢！請查看當前資料夾。');


% =========================================================================
% 5. 繪製推論曲面圖
% =========================================================================

% 建立一個較大的視窗來容納 4 個子圖(時間)
figure('Name', '去模糊化方法比較 - 加熱時間 (Time_z)', 'Position', [100, 100, 1000, 800]);

% --- 1. COG (重心法) ---
subplot(2, 2, 1);
surf(X, Y, reshape(out_cog(:, 1), size(X)));
title('COG (重心法)');
xlabel('Temp x (°C)'); ylabel('Weight y (kg)'); zlabel('Time (min)');
view(45, 30); % 設定視角

% --- 2. MOM (最大平均法) ---
subplot(2, 2, 2);
surf(X, Y, reshape(out_mom(:, 1), size(X)));
title('MOM (最大平均法)');
xlabel('Temp x (°C)'); ylabel('Weight y (kg)'); zlabel('Time (min)');
view(45, 30);

% --- 3. mMOM (修改型最大平均法) ---
subplot(2, 2, 3);
surf(X, Y, reshape(out_mmom(:, 1), size(X)));
title('mMOM (修改型最大平均法)');
xlabel('Temp x (°C)'); ylabel('Weight y (kg)'); zlabel('Time (min)');
view(45, 30);

% --- 4. CA (中心平均法) ---
subplot(2, 2, 4);
surf(X, Y, reshape(out_ca(:, 1), size(X)));
title('CA (中心平均法)');
xlabel('Temp x (°C)'); ylabel('Weight y (kg)'); zlabel('Time (min)');
view(45, 30);

all_axes = findobj(gcf, 'Type', 'axes');

% 一次性為所有抓到的座標軸設定 X, Y, Z 刻度
set(all_axes, 'XTick', -4:1:4);
set(all_axes, 'YTick', 0:0.5:2);
set(all_axes, 'ZTick', 0:3:15);

sgtitle('四種去模糊化方法之控制曲面比較 (輸出：加熱時間)');


% 建立一個較大的視窗來容納 4 個子圖(功率)
figure('Name', '去模糊化方法比較 - Power (Watt)', 'Position', [100, 100, 1000, 800]);

% --- 1. COG (重心法) ---
subplot(2, 2, 1);
surf(X, Y, reshape(out_cog(:, 2), size(X)));
title('COG (重心法)');
xlabel('Temp x (°C)'); ylabel('Weight y (kg)'); zlabel('Power (W)');
view(45, 30); % 設定視角

% --- 2. MOM (最大平均法) ---
subplot(2, 2, 2);
surf(X, Y, reshape(out_mom(:, 2), size(X)));
title('MOM (最大平均法)');
xlabel('Temp x (°C)'); ylabel('Weight y (kg)'); zlabel('Power (W)');
view(45, 30);

% --- 3. mMOM (修改型最大平均法) ---
subplot(2, 2, 3);
surf(X, Y, reshape(out_mmom(:, 2), size(X)));
title('mMOM (修改型最大平均法)');
xlabel('Temp x (°C)'); ylabel('Weight y (kg)'); zlabel('Power (W)');
view(45, 30);

% --- 4. CA (中心平均法) ---
subplot(2, 2, 4);
surf(X, Y, reshape(out_ca(:, 2), size(X)));
title('CA (中心平均法)');
xlabel('Temp x (°C)'); ylabel('Weight y (kg)'); zlabel('Power (W)');
view(45, 30);

all_axes = findobj(gcf, 'Type', 'axes');

% 一次性為所有抓到的座標軸設定 X, Y, Z 刻度
set(all_axes, 'XTick', -4:1:4);
set(all_axes, 'YTick', 0:0.5:2);
set(all_axes, 'ZTick', 600:100:1200);

sgtitle('四種去模糊化方法之控制曲面比較 (輸出：微波功率)');