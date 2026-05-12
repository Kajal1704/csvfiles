import matplotlib
matplotlib.use("Agg")

from fastapi import FastAPI, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
import pandas as pd
import numpy as np
import joblib
import shap
import matplotlib.pyplot as plt
import base64
from io import BytesIO
import io
import os
from xgboost import XGBClassifier

app = FastAPI()

# ============================
# ENABLE CORS
# ============================
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ============================
# LOAD MODEL
# ============================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

MODEL_PATH = os.path.join(BASE_DIR, "xgboost_disease_model.json")
ENCODER_PATH = os.path.join(BASE_DIR, "label_encoder.pkl")
FEATURE_PATH = os.path.join(BASE_DIR, "feature_columns.pkl")

model = XGBClassifier()
model.load_model(MODEL_PATH)

label_encoder = joblib.load(ENCODER_PATH)
feature_columns = joblib.load(FEATURE_PATH)

print("Model loaded successfully")


# ============================
# HOME ROUTE
# ============================
@app.get("/")
def home():
    return {"message": "FastAPI backend running"}


# ============================
# FEATURE EXTRACTION (PIPELINE ALIGNED)
# ============================
def extract_features(df):

    # -------- GNN SIMULATION (Mutation Detection) --------
    mutation_detected = int(df["Mutation"].max() > 0.3)

    # -------- Pathogenic Detection --------
    #pathogenic_score = (
     #   0.6 * df["Conservation_Score"].mean() +
     #   0.4 * df["Expression"].mean()
   # )

   

    #clinical_significance = 1 if pathogenic_score > 0.7 else 0
    pathogenic_score = (
    0.6 * df["Conservation_Score"].mean() +
    0.4 * df["Expression"].mean()
)

    # Pathogenicity only if mutation exists
    if mutation_detected == 1:
        clinical_significance = 1 if pathogenic_score > 0.7 else 0
    else:
        clinical_significance = 0

    # -------- Risk Score --------
    #risk_score = float(df["Expression"].mean())
    risk_score = (
    0.4 * df["Expression"].mean()
    + 0.3 * df["Mutation"].mean()
    + 0.3 * df["Conservation_Score"].mean())

    return [risk_score, mutation_detected, clinical_significance]


# ============================
# PREDICTION API
# ============================
@app.post("/predict")
async def predict(file: UploadFile = File(...)):

    try:
        contents = await file.read()
        input_df = pd.read_csv(io.BytesIO(contents))

        # Extract features
        features = extract_features(input_df)
        X = pd.DataFrame([features], columns=feature_columns)

        # -------- PIPELINE LOGIC --------

        # Extract features clearly
        risk_score = float(features[0])
        mutation_detected = int(features[1])
        clinical_significance = int(features[2])

        # 1️⃣ No Mutation
        if mutation_detected == 0:
            disease = "No Mutation Detected"
            risk_level = "Low Risk"

        # 2️⃣ Benign Mutation
        elif mutation_detected==1 and clinical_significance == 0:
            disease = "Benign Mutation (No Disease)"
            risk_level = "Low Risk"

        # 3️⃣ Pathogenic → Disease Prediction
        else:
            genes = list(input_df["Gene_Symbol"].unique())

            
            cf_count = genes.count("CFTR")

            scd_count = genes.count("HBB") + genes.count("HBD")

            hd_count = genes.count("HTT")
        
            if risk_score < 0.8:
                risk_level = "Low Risk"
            elif risk_score < 1.5:
                risk_level = "Medium Risk"
            else:
                risk_level = "High Risk"

# Predict dominant disease
            if hd_count > scd_count and hd_count > cf_count:
                 disease = "Huntington's Disease (HD)"
                

            elif scd_count > hd_count and scd_count > cf_count:
                disease = "Sickle Cell Disease (SCD)"
                

            elif cf_count > hd_count and cf_count > scd_count:
                disease = "Cystic Fibrosis (CF)"
                
            else:
                prediction = int(model.predict(X)[0])
                disease = str(label_encoder.inverse_transform([prediction])[0])
                

        # Final response
        return {
            "prediction": disease,
            "risk_level": risk_level,
            "features_used": {
            "Risk_Score": risk_score,
            "Mutation_Status": mutation_detected,
            "Clinical_Significance": clinical_significance
        }
    }

    except Exception as e:
        return {"error": str(e)}


# ============================
# GLOBAL SHAP
# ============================
@app.get("/explain")
def explain():

    data_path = os.path.join(BASE_DIR, "rf_training_dataset.csv")
    df = pd.read_csv(data_path)
    X = df[feature_columns]

    explainer = shap.Explainer(model, X)
    shap_values = explainer(X)

    shap_values_mean = np.abs(shap_values.values).mean(axis=2)

    # Feature Importance Plot
    plt.figure()
    shap.summary_plot(
        shap_values_mean,
        X,
        plot_type="bar",
        show=False
    )

    buf1 = BytesIO()
    plt.savefig(buf1, format="png", bbox_inches="tight")
    plt.close()
    buf1.seek(0)

    global_plot = base64.b64encode(buf1.read()).decode()

    # Beeswarm Plot
    plt.figure()
    shap.summary_plot(
        shap_values_mean,
        X,
        show=False
    )

    buf2 = BytesIO()
    plt.savefig(buf2, format="png", bbox_inches="tight")
    plt.close()
    buf2.seek(0)

    summary_plot = base64.b64encode(buf2.read()).decode()

    return {
        "global_plot": global_plot,
        "summary_plot": summary_plot
    }


# ============================
# LOCAL SHAP
# ============================
@app.post("/explain-local")
async def explain_local(file: UploadFile = File(...)):

    try:
        contents = await file.read()
        df = pd.read_csv(io.BytesIO(contents))

        # Extract same features used in prediction
        features = extract_features(df)
        X = pd.DataFrame([features], columns=feature_columns)

        explainer = shap.Explainer(model)
        shap_values = explainer(X)

        prediction = int(model.predict(X)[0])

        shap_vals = shap_values.values[0][prediction]
        base_val = explainer.expected_value[prediction]

        force_plot = shap.force_plot(
            base_val,
            shap_vals,
            X.iloc[0],
            feature_names=feature_columns
        )

        html_buffer = io.StringIO()
        shap.save_html(html_buffer, force_plot)

        return {"html": html_buffer.getvalue()}

    except Exception as e:
        return {"error": str(e)}


# ============================
# MODEL COMPARISON
# ============================
@app.get("/model-comparison")
def model_comparison():

    try:
        import seaborn as sns
        from sklearn.metrics import confusion_matrix
        from sklearn.ensemble import RandomForestClassifier

        data_path = os.path.join(BASE_DIR, "rf_training_dataset.csv")
        df = pd.read_csv(data_path)
        X = df[feature_columns]
        y = label_encoder.transform(df["disease"])

        # XGBoost
        xgb_preds = model.predict(X)
        xgb_cm = confusion_matrix(y, xgb_preds)
        xgb_acc = (xgb_preds == y).mean()

        plt.figure()
        sns.heatmap(xgb_cm, annot=True, fmt="d", cmap="Blues")
        plt.title("XGBoost Confusion Matrix")

        buf = BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)

        xgb_img = base64.b64encode(buf.read()).decode("utf-8")
        plt.close()

        # Random Forest
        rf = RandomForestClassifier(random_state=42)
        rf.fit(X, y)

        rf_preds = rf.predict(X)
        rf_cm = confusion_matrix(y, rf_preds)
        rf_acc = (rf_preds == y).mean()

        plt.figure()
        sns.heatmap(rf_cm, annot=True, fmt="d", cmap="Greens")
        plt.title("Random Forest Confusion Matrix")

        buf = BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)

        rf_img = base64.b64encode(buf.read()).decode("utf-8")
        plt.close()

        return {
            "rf_accuracy": round(rf_acc * 100, 2),
            "xgb_accuracy": round(xgb_acc * 100, 2),
            "rf_confusion_matrix": rf_img,
            "xgb_confusion_matrix": xgb_img
        }

    except Exception as e:
        return {"error": str(e)}