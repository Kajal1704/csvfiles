// import { useState } from "react";

// export default function Predict() {

//   const [file, setFile] = useState(null);
//   const [prediction, setPrediction] = useState(null);

//   const handlePredict = async () => {
//     if (!file) {
//       alert("Please upload a CSV file");
//       return;
//     }

//     const formData = new FormData();
//     formData.append("file", file);

//     try {
//       const response = await fetch("http://127.0.0.1:8000/predict", {
//         method: "POST",
//         body: formData,
//       });

//       const data = await response.json();
//       setPrediction(data);

//     } catch (error) {
//       console.error("Prediction error:", error);
//       alert("Error connecting to backend");
//     }
//   };

//   const handleReset = () => {
//     setFile(null);
//     setPrediction(null);
//   };

//   return (
//     <div className="min-h-screen bg-[#f5f7f9]">

//       <div className="px-20 py-12">

//         {/* Title */}
//         <h2 className="text-3xl font-bold text-gray-800">
//           Disease Prediction
//         </h2>

//         <p className="text-gray-500 mt-2">
//           Upload a genetic CSV file to predict the associated disease.
//         </p>

//         {/* Two Column Layout */}
//         <div className="grid grid-cols-2 gap-8 mt-10">

//           {/* LEFT CARD */}
//           <div className="bg-white rounded-xl shadow-sm p-8">

//             <h3 className="text-lg font-semibold text-gray-700 mb-6">
//               Input Features
//             </h3>

//             {/* CSV Upload */}
//             <div className="mb-6">
//               <label className="block text-sm font-medium mb-2">
//                 Upload the CSV File
//               </label>

//               <input
//                 type="file"
//                 accept=".csv"
//                 onChange={(e) => setFile(e.target.files[0])}
//                 className="w-full border rounded p-2"
//               />
//             </div>

//             {/* Buttons */}
//             <div className="flex gap-4">
//               <button
//                 onClick={handlePredict}
//                 className="w-full bg-[#6c8ea3] text-white py-2 rounded-lg hover:opacity-90 transition"
//               >
//                 Predict Disease
//               </button>

//               <button
//                 onClick={handleReset}
//                 className="border border-gray-300 px-6 py-2 rounded-lg text-gray-600 hover:bg-gray-100 transition"
//               >
//                 Reset
//               </button>
//             </div>

//           </div>

//           {/* RIGHT CARD */}
//           <div className="bg-white rounded-xl shadow-sm p-8 flex flex-col justify-center items-center text-center">

//             {prediction && prediction.predictions ? (
//               <>
//                 <div className="text-5xl mb-4">🧬</div>

//                 <h4 className="text-lg font-semibold text-gray-700">
//                   Prediction Results
//                 </h4>

//                 <div className="mt-6 space-y-2">
//                   {prediction.predictions.map((item, index) => (
//   <div key={index} className="mb-4">
//     <div className="text-lg font-semibold text-gray-800">
//       Disease: {item.predicted_disease}
//     </div>

//     <div className="text-sm text-gray-500">
//       Confidence: {item.confidence}%
//     </div>
//   </div>
// ))}
//                 </div>
//               </>
//             ) : (
//               <>
//                 <div className="text-5xl text-gray-300 mb-4">🧪</div>

//                 <h4 className="text-lg font-semibold text-gray-700">
//                   No prediction yet
//                 </h4>

//                 <p className="text-gray-500 mt-2">
//                   Upload a CSV file and click "Predict Disease"
//                 </p>
//               </>
//             )}

//           </div>

//         </div>

//       </div>

//       {/* Footer */}
//       <footer className="text-center text-sm text-gray-500 py-6 border-t mt-12">
//         © 2026 Genetic Mutation Detection & Disease Prediction System — Academic Research Project
//       </footer>

//     </div>
//   );
// }


// import { useState } from "react";

// export default function Predict() {

//   const [file, setFile] = useState(null);
//   const [prediction, setPrediction] = useState(null);
//   const [loading, setLoading] = useState(false);

//   const handlePredict = async () => {

//     if (!file) {
//       alert("Please upload a CSV file");
//       return;
//     }

//     const formData = new FormData();
//     formData.append("file", file);

//     try {

//       setLoading(true);

//       const response = await fetch("http://127.0.0.1:8000/predict", {
//         method: "POST",
//         body: formData,
//       });

//       const data = await response.json();

//       console.log("Backend Response:", data);

//       if (data.error) {
//         alert(data.error);
//         return;
//       }

//       setPrediction(data);

//     } catch (error) {

//       console.error("Prediction error:", error);
//       alert("Error connecting to backend");

//     } finally {
//       setLoading(false);
//     }
//   };

//   const handleReset = () => {
//     setFile(null);
//     setPrediction(null);
//   };

//   return (
//     <div className="min-h-screen bg-[#f5f7f9]">

//       <div className="px-20 py-12">

//         <h2 className="text-3xl font-bold text-gray-800">
//           Disease Prediction
//         </h2>

//         <p className="text-gray-500 mt-2">
//           Upload a genetic CSV file to predict the associated disease.
//         </p>

//         <div className="grid grid-cols-2 gap-8 mt-10">

//           {/* LEFT CARD */}
//           <div className="bg-white rounded-xl shadow-sm p-8">

//             <h3 className="text-lg font-semibold text-gray-700 mb-6">
//               Input Features
//             </h3>

//             <div className="mb-6">

//               <label className="block text-sm font-medium mb-2">
//                 Upload the CSV File
//               </label>

//               <input
//                 type="file"
//                 accept=".csv"
//                 onChange={(e) => setFile(e.target.files[0])}
//                 className="w-full border rounded p-2"
//               />

//             </div>

//             <div className="flex gap-4">

//               <button
//                 onClick={handlePredict}
//                 className="w-full bg-[#6c8ea3] text-white py-2 rounded-lg hover:opacity-90 transition"
//               >
//                 {loading ? "Predicting..." : "Predict Disease"}
//               </button>

//               <button
//                 onClick={handleReset}
//                 className="border border-gray-300 px-6 py-2 rounded-lg text-gray-600 hover:bg-gray-100 transition"
//               >
//                 Reset
//               </button>

//             </div>

//           </div>


//           {/* RIGHT CARD */}
//           <div className="bg-white rounded-xl shadow-sm p-8 flex flex-col justify-center items-center text-center">

//             {prediction ? (

//               <>
//                 <div className="text-5xl mb-4">🧬</div>

//                 <h4 className="text-lg font-semibold text-gray-700">
//                   Prediction Results
//                 </h4>

//                 <div className="mt-6 space-y-2">

//                   <div className="text-xl font-bold text-green-700">
//                     {prediction.prediction}
//                   </div>

//                   <div className="text-sm text-gray-500">
//                     Risk Score: {prediction.features_used.Risk_Score}
//                   </div>

//                   <div className="text-sm text-gray-500">
//                     Mutation Status: {prediction.features_used.Mutation_Status}
//                   </div>

//                   <div className="text-sm text-gray-500">
//                     Clinical Significance: {prediction.features_used.Clinical_Significance}
//                   </div>

//                 </div>
//               </>

//             ) : (

//               <>
//                 <div className="text-5xl text-gray-300 mb-4">🧪</div>

//                 <h4 className="text-lg font-semibold text-gray-700">
//                   No prediction yet
//                 </h4>

//                 <p className="text-gray-500 mt-2">
//                   Upload a CSV file and click "Predict Disease"
//                 </p>
//               </>

//             )}

//           </div>

//         </div>

//       </div>

//       <footer className="text-center text-sm text-gray-500 py-6 border-t mt-12">
//         © 2026 Genetic Mutation Detection & Disease Prediction System — Academic Research Project
//       </footer>

//     </div>
//   );
// }


import { useState } from "react";

export default function Predict() {

  const [file, setFile] = useState(null);
  const [prediction, setPrediction] = useState(null);
  const [loading, setLoading] = useState(false);

  // Map numeric prediction → disease name
  const diseaseMap = {
    0: "Cystic Fibrosis (CF)",
    1: "Sickle Cell Disease (SCD)",
    2: "Huntington’s Disease (HD)"
  };

  const handlePredict = async () => {

    if (!file) {
      alert("Please upload a CSV file");
      return;
    }

    const formData = new FormData();
    formData.append("file", file);

    try {

      setLoading(true);

      const response = await fetch("http://127.0.0.1:8000/predict", {
        method: "POST",
        body: formData,
      });

      const data = await response.json();

      console.log("Backend Response:", data);

      if (data.error) {
        alert(data.error);
        return;
      }

      setPrediction(data);

    } catch (error) {

      console.error("Prediction error:", error);
      alert("Error connecting to backend");

    } finally {
      setLoading(false);
    }
  };

  const handleReset = () => {
    setFile(null);
    setPrediction(null);
  };

  return (
    <div className="min-h-screen bg-[#f5f7f9]">

      <div className="px-20 py-12">

        {/* Title */}
        <h2 className="text-3xl font-bold text-gray-800">
          Disease Prediction
        </h2>

        <p className="text-gray-500 mt-2">
          Upload a genetic CSV file to predict the associated disease.
        </p>

        <div className="grid grid-cols-2 gap-8 mt-10">

          {/* LEFT CARD */}
          <div className="bg-white rounded-xl shadow-sm p-8">

            <h3 className="text-lg font-semibold text-gray-700 mb-6">
              Input Features
            </h3>

            <div className="mb-6">

              <label className="block text-sm font-medium mb-2">
                Upload the CSV File
              </label>

              <input
                type="file"
                accept=".csv"
                onChange={(e) => setFile(e.target.files[0])}
                className="w-full border rounded p-2"
              />

            </div>

            <div className="flex gap-4">

              <button
                onClick={handlePredict}
                className="w-full bg-[#6c8ea3] text-white py-2 rounded-lg hover:opacity-90 transition"
              >
                {loading ? "Predicting..." : "Predict Disease"}
              </button>

              <button
                onClick={handleReset}
                className="border border-gray-300 px-6 py-2 rounded-lg text-gray-600 hover:bg-gray-100 transition"
              >
                Reset
              </button>

            </div>

          </div>


          {/* RIGHT CARD */}
          <div className="bg-white rounded-xl shadow-sm p-8 flex flex-col justify-center items-center text-center">

            {prediction ? (

              <>
                <div className="text-5xl mb-4">🧬</div>

                <h4 className="text-lg font-semibold text-gray-700">
                  Prediction Results
                </h4>

                <div className="mt-6 space-y-2">

                  <div className="text-xl font-bold text-green-700">

                    {diseaseMap[prediction.prediction] || prediction.prediction}

                  </div>

                  <div className="text-sm text-gray-500">
                    Risk Score: {prediction.features_used.Risk_Score}
                  </div>

                  <div className="text-sm text-gray-500">
                    Mutation Status: {prediction.features_used.Mutation_Status}
                  </div>

                  <div className="text-sm text-gray-500">
                    Clinical Significance: {prediction.features_used.Clinical_Significance}
                  </div>

                </div>

              </>

            ) : (

              <>
                <div className="text-5xl text-gray-300 mb-4">🧪</div>

                <h4 className="text-lg font-semibold text-gray-700">
                  No prediction yet
                </h4>

                <p className="text-gray-500 mt-2">
                  Upload a CSV file and click "Predict Disease"
                </p>
              </>

            )}

          </div>

        </div>

      </div>

      <footer className="text-center text-sm text-gray-500 py-6 border-t mt-12">
        © 2026 Genetic Mutation Detection & Disease Prediction System — Academic Research Project
      </footer>

    </div>
  );
}