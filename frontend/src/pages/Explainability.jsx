import React, { useEffect, useState } from "react";

const Explainability = () => {

  const [globalPlots, setGlobalPlots] = useState(null);
  const [shapHtml, setShapHtml] = useState(null);
  const [file, setFile] = useState(null);
  const [loading, setLoading] = useState(true);

  // ================= GLOBAL SHAP =================
  useEffect(() => {

    fetch("http://127.0.0.1:8000/explain")
      .then((res) => res.json())
      .then((data) => {
        setGlobalPlots(data);
        setLoading(false);
      })
      .catch((err) => {
        console.error("Global SHAP error:", err);
        setLoading(false);
      });

  }, []);

  // ================= LOCAL SHAP =================
  const handleLocalExplain = async () => {

    if (!file) {
      alert("Please upload CSV first");
      return;
    }

    const formData = new FormData();
    formData.append("file", file);

    try {

      const response = await fetch(
        "http://127.0.0.1:8000/explain-local",
        {
          method: "POST",
          body: formData
        }
      );

      const data = await response.json();

      if (data.error) {
        alert(data.error);
        return;
      }

      setShapHtml(data.html);

    } catch (error) {
      console.error("Local SHAP error:", error);
    }

  };

  return (

    <div className="min-h-screen bg-[#f5f7f9] px-20 py-12">

      {/* HEADER */}
      <div className="mb-10 flex justify-between items-start">

        <div>
          <h1 className="text-3xl font-bold text-gray-800">
            Model Explainability
          </h1>

          <p className="text-gray-500 mt-2">
            This section helps understand how the system makes predictions from genetic data.
          </p>
        </div>

        <div className="bg-green-100 text-green-700 px-4 py-1 rounded-full text-sm font-medium">
          Model: XGBoost | SHAP
        </div>

      </div>


      {/* GLOBAL EXPLAINABILITY */}
      {loading ? (

        <div className="text-center text-gray-500">
          Loading visualizations...
        </div>

      ) : globalPlots ? (

        <>
          <h2 className="text-2xl font-semibold mb-4">
            Overall Model Behavior
          </h2>

          {/* SIMPLE EXPLANATION */}
          <div className="bg-green-50 border-l-4 border-green-500 p-5 rounded-lg mb-8">

            <p className="text-gray-700 text-base">
              These graphs show how the system makes decisions across many patients.
            </p>

            <p className="text-gray-700 text-base leading-relaxed mt-2">
              The model mainly relies on <strong>Risk Score</strong>, which represents
              how severe or important a mutation is.
            </p>

            <p className="text-gray-700 text-base leading-relaxed mt-2">
              It also considers whether a mutation is harmful (<strong>Clinical Significance</strong>)
              and whether a mutation is present (<strong>Mutation Status</strong>).
            </p>

            <p className="text-gray-700 text-base leading-relaxed mt-2">
              Together, these help the system decide whether a disease is likely or not.
            </p>

          </div>

          {/* SHAP PLOTS */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-8 mb-10">

            <div className="bg-white rounded-xl shadow-sm p-6">
              <h3 className="text-lg font-semibold mb-4">
                Important Factors in Prediction
              </h3>
              <img
                src={`data:image/png;base64,${globalPlots.global_plot}`}
                alt="Feature Importance"
                className="rounded-xl w-full"
              />
            </div>

            <div className="bg-white rounded-xl shadow-sm p-6">
              <h3 className="text-lg font-semibold mb-4">
                Feature Impact Distribution
              </h3>
              <img
                src={`data:image/png;base64,${globalPlots.summary_plot}`}
                alt="Summary Plot"
                className="rounded-xl w-full"
              />
            </div>

          </div>

          {/* SIMPLE SUMMARY CARDS */}
          <div className="grid grid-cols-3 gap-6 mb-12">

            <div className="bg-[#f1f5f8] p-4 rounded-xl text-center">
              <p className="text-xs text-gray-500">Main Factor</p>
              <p className="font-semibold text-lg">Risk Score</p>
            </div>

            <div className="bg-[#f1f5f8] p-4 rounded-xl text-center">
              <p className="text-xs text-gray-500">Helps Identify Harm</p>
              <p className="font-semibold text-lg">Clinical Significance</p>
            </div>

            <div className="bg-[#f1f5f8] p-4 rounded-xl text-center">
              <p className="text-xs text-gray-500">Indicates Presence</p>
              <p className="font-semibold text-lg">Mutation Status</p>
            </div>

          </div>

        </>

      ) : (

        <div className="text-red-500">
          Failed to load visualizations.
        </div>

      )}


      {/* LOCAL EXPLAINABILITY */}
      <h2 className="text-2xl font-semibold mb-4">
        Individual Prediction Explanation
      </h2>

      {/* SIMPLE EXPLANATION */}
      <div className="bg-green-50 border-l-4 border-green-500 p-5 rounded-lg mb-6">

        <p className="text-gray-700 text-base">
          This section explains the prediction for a single uploaded patient file.
        </p>

        <p className="text-gray-700 text-base leading-relaxed mt-2">
          Each factor either increases or decreases the chance of disease.
        </p>

        <div className="mt-3 text-base">
          <span className="text-red-500 font-semibold">
            Red → pushes toward disease
          </span><br/>
          <span className="text-blue-500 font-semibold">
            Blue → pushes toward normal (no disease)
          </span>
        </div>

        <p className="text-gray-600 text-base leading-relaxed mt-2">
          Larger bars indicate stronger influence on the final decision.
        </p>

      </div>


      <div className="bg-white rounded-xl shadow-sm p-8">

        <h3 className="text-lg font-semibold mb-4">
          SHAP Force Plot
        </h3>

        {/* Upload */}
        <input
          type="file"
          accept=".csv"
          onChange={(e)=>setFile(e.target.files[0])}
          className="mb-4"
        />

        <button
          onClick={handleLocalExplain}
          className="bg-[#6c8ea3] text-white px-4 py-2 rounded-lg mb-6 hover:opacity-90"
        >
          Generate Explanation
        </button>

        {/* SHAP Plot */}
        {shapHtml ? (

          <iframe
            title="SHAP Force Plot"
            srcDoc={shapHtml}
            style={{
              width: "100%",
              height: "400px",
              border: "none"
            }}
          />

        ) : (

          <div className="bg-gray-100 rounded-xl h-64 flex items-center justify-center text-gray-400">
            Upload a file to view explanation
          </div>

        )}

        {/* INTERPRETATION */}
        <div className="bg-yellow-50 border border-yellow-300 p-5 rounded-lg mt-6">
          <h3 className="text-md font-semibold mb-2">
            Clinical Interpretation
          </h3>

          <p className="text-base text-gray-700">
            In this prediction, the system mainly looks at how severe the mutation is
            (Risk Score), whether it is harmful (Clinical Significance),
            and whether it exists (Mutation Status).
          </p>

          <p className="text-base text-gray-700 mt-2">
            These combined factors help determine whether a disease is likely.
          </p>

        </div>

      </div>

    </div>
  );
};

export default Explainability;