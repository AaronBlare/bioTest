import pandas as pd
import numpy as np
import re
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error
from scipy.stats import pearsonr
import lightgbm as lgb
import warnings
warnings.filterwarnings('ignore')

# ------------------------------
# Функция для очистки имен признаков
# ------------------------------
def clean_feature_name(name):
    """
    Заменяет все символы, кроме букв, цифр и подчеркивания, на подчеркивание.
    Также удаляет начальные/конечные подчеркивания, если они появились.
    """
    # Заменяем все не-буквенно-цифровые и не-подчеркивание на подчеркивание
    cleaned = re.sub(r'[^\w]', '_', str(name))
    # Убираем множественные подчеркивания подряд
    cleaned = re.sub(r'_+', '_', cleaned)
    # Убираем подчеркивания в начале и конце
    cleaned = cleaned.strip('_')
    # Если имя стало пустым, даем дефолтное имя
    if not cleaned:
        cleaned = 'feature'
    return cleaned

# ------------------------------
# 1. Загрузка данных
# ------------------------------
# Читаем список признаков из файла feats.xlsx
feats_df = pd.read_excel('D:/Yandex.Disk/bbd/sex_hormones/ДНКОМ/feats_sets/4_1/feats.xlsx')
features_raw = feats_df['Features'].dropna().tolist()

# Очищаем имена признаков
features_clean = [clean_feature_name(f) for f in features_raw]
# Сохраняем маппинг для возможного переименования в data.xlsx
feature_mapping = dict(zip(features_raw, features_clean))

print(f"Загружено признаков из feats.xlsx: {len(features_raw)}")
print("Пример очистки:")
for orig, clean in list(feature_mapping.items())[:3]:
    print(f"  '{orig}' -> '{clean}'")

# Читаем основные данные
data = pd.read_excel('D:/Yandex.Disk/bbd/sex_hormones/ДНКОМ/feats_sets/4_1/data.xlsx')
print(f"\nРазмер data.xlsx: {data.shape}")

# Проверяем наличие целевой переменной
if 'Возраст' not in data.columns:
    raise ValueError("В data.xlsx отсутствует столбец 'Возраст'")

# Переименовываем столбцы в data.xlsx согласно очищенным именам
# Сначала находим, какие исходные признаки присутствуют в data
existing_original = [f for f in features_raw if f in data.columns]
# Для них применяем переименование
rename_dict = {orig: feature_mapping[orig] for orig in existing_original}
data.rename(columns=rename_dict, inplace=True)

# Теперь список используемых признаков - это очищенные имена
features_to_use = [feature_mapping[orig] for orig in existing_original]

print(f"Используется признаков (после очистки): {len(features_to_use)}")
if len(features_to_use) == 0:
    raise ValueError("Нет ни одного признака из feats.xlsx в data.xlsx после переименования!")

# ------------------------------
# 2. Подготовка X и y
# ------------------------------
X = data[features_to_use].copy()
y = data['Возраст'].copy()

# Проверка на пропуски
if X.isnull().any().any():
    print("\nОбнаружены пропуски в признаках. Заполняем медианой...")
    X = X.fillna(X.median())
if y.isnull().any():
    print("Обнаружены пропуски в целевой переменной. Удаляем такие строки.")
    mask = ~y.isnull()
    X = X[mask]
    y = y[mask]

print(f"\nФинальный размер выборки: {X.shape[0]} образцов, {X.shape[1]} признаков")

# ------------------------------
# 3. Разделение на train / validation / test (60/20/20)
# ------------------------------
X_train, X_temp, y_train, y_temp = train_test_split(
    X, y, test_size=0.4, random_state=42, shuffle=True
)
X_val, X_test, y_val, y_test = train_test_split(
    X_temp, y_temp, test_size=0.5, random_state=42, shuffle=True
)

print(f"Размер обучающей выборки: {len(X_train)}")
print(f"Размер валидационной выборки: {len(X_val)}")
print(f"Размер тестовой выборки: {len(X_test)}")

# ------------------------------
# 4. Обучение модели LightGBM
# ------------------------------
params = {
    'objective': 'regression',
    'metric': 'mae',
    'boosting_type': 'gbdt',
    'num_leaves': 31,
    'learning_rate': 0.05,
    'feature_fraction': 0.8,
    'bagging_fraction': 0.8,
    'bagging_freq': 5,
    'verbose': -1,
    'random_state': 42
}

train_data = lgb.Dataset(X_train, label=y_train)
val_data = lgb.Dataset(X_val, label=y_val, reference=train_data)

model = lgb.train(
    params,
    train_data,
    valid_sets=[val_data],
    num_boost_round=1000,
    callbacks=[lgb.early_stopping(stopping_rounds=50), lgb.log_evaluation(100)]
)

# ------------------------------
# 5. Оценка метрик
# ------------------------------
def evaluate(y_true, y_pred, set_name=""):
    mae = mean_absolute_error(y_true, y_pred)
    corr, p_value = pearsonr(y_true, y_pred)
    print(f"\n{set_name} результаты:")
    print(f"  MAE = {mae:.3f}")
    print(f"  Корреляция Пирсона = {corr:.4f} (p-value = {p_value:.2e})")
    return mae, corr

y_val_pred = model.predict(X_val)
evaluate(y_val, y_val_pred, "Валидационная")

y_test_pred = model.predict(X_test)
evaluate(y_test, y_test_pred, "Тестовая")

# ------------------------------
# 6. Важность признаков
# ------------------------------
importance = model.feature_importance(importance_type='gain')
feature_importance_df = pd.DataFrame({
    'feature': features_to_use,
    'importance': importance
}).sort_values('importance', ascending=False)

print("\nТоп-4 наиболее важных признаков по gain:")
print(feature_importance_df.head(4).to_string(index=False))

# Сохраняем модель (опционально)
# model.save_model('age_regression_model.txt')
print("\nМодель успешно обучена и оценена.")