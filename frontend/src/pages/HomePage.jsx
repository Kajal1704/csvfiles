// import React from "react";
// import { Link, useLocation } from "react-router-dom";
// import {
//   Home,
//   FlaskConical,
//   Brain,
//   BarChart3,
//   ShieldCheck,
//   ArrowRight,
// } from "lucide-react";

// export default function HomePage() {
//   const location = useLocation();

//   return (
//     <div className="min-h-screen bg-gray-50">

//       {/* ================= NAVBAR ================= */}
//       <nav className="bg-white border-b border-gray-200 px-12 py-4 flex justify-between items-center">
        
//         <div className="flex items-center gap-2 text-teal-700 font-semibold text-lg">
//           🧬 <span>GeneMutate AI</span>
//         </div>

//         <div className="flex items-center gap-8 text-sm font-medium text-gray-600">
          
//           {/* HOME */}
//           <Link
//             to="/"
//             className={`flex items-center gap-2 px-4 py-2 rounded-lg transition ${
//               location.pathname === "/"
//                 ? "bg-teal-600 text-white"
//                 : "hover:text-teal-600"
//             }`}
//           >
//             <Home size={16} />
//             Home
//           </Link>

//           {/* PREDICT */}
//           <Link
//             to="/predict"
//             className={`transition ${
//               location.pathname === "/predict"
//                 ? "text-teal-600 font-semibold"
//                 : "hover:text-teal-600"
//             }`}
//           >
//             Predict
//           </Link>

//           <button className="hover:text-teal-600 transition">
//             Explainability
//           </button>

//           <button className="hover:text-teal-600 transition">
//             Model Comparison
//           </button>

//           <button className="hover:text-teal-600 transition">
//             About
//           </button>
//         </div>
//       </nav>

//       {/* ================= HERO SECTION ================= */}
//       <section className="bg-gradient-to-r from-teal-700 to-emerald-600 text-white text-center py-32 px-6">
//         <div className="max-w-4xl mx-auto">

//           <div className="inline-flex items-center gap-2 bg-white/20 px-6 py-2 rounded-full text-sm backdrop-blur-md mb-8">
//             🔬 Machine Learning–Powered Genomic Analysis
//           </div>

//           <h1 className="text-5xl font-bold leading-tight mb-8">
//             Genetic Mutation Detection &
//             <br />
//             Disease Prediction System
//           </h1>

//           <p className="text-lg text-gray-100 leading-relaxed mb-10 max-w-3xl mx-auto">
//             Leveraging XGBoost and Explainable AI (SHAP) to detect genetic
//             mutations and classify diseases — Sickle Cell Disease,
//             Cystic Fibrosis, and Huntington's Disease — with
//             transparent, interpretable predictions.
//           </p>

//           <div className="flex justify-center gap-6">
            
//             {/* START PREDICTION BUTTON */}
//             <Link
//               to="/predict"
//               className="flex items-center gap-2 bg-emerald-400 text-white px-8 py-3 rounded-xl font-semibold hover:opacity-90 transition"
//             >
//               <FlaskConical size={18} />
//               Start Prediction
//             </Link>

//             <button className="flex items-center gap-2 border border-white px-8 py-3 rounded-xl font-semibold hover:bg-white hover:text-teal-700 transition">
//               Learn More
//               <ArrowRight size={16} />
//             </button>

//           </div>
//         </div>
//       </section>

//       {/* ================= HOW IT WORKS ================= */}
//       <section className="bg-gray-100 py-24 px-12">
//         <div className="text-center mb-16">
//           <h2 className="text-3xl font-bold text-gray-800 mb-3">
//             How It Works
//           </h2>
//           <p className="text-gray-500">
//             A step-by-step overview of the prediction pipeline
//           </p>
//         </div>

//         <div className="grid md:grid-cols-4 gap-10 max-w-7xl mx-auto">

//           <StepCard
//             icon={<FlaskConical size={26} />}
//             step="STEP 1"
//             title="Input Genetic Data"
//             desc="Enter mutation features: risk score, mutation status, and clinical significance."
//           />

//           <StepCard
//             icon={<Brain size={26} />}
//             step="STEP 2"
//             title="ML Prediction"
//             desc="XGBoost model classifies the mutation into SCD, CF, or Huntington's Disease."
//           />

//           <StepCard
//             icon={<BarChart3 size={26} />}
//             step="STEP 3"
//             title="SHAP Explanation"
//             desc="Explainable AI reveals which features drove the prediction."
//           />

//           <StepCard
//             icon={<ShieldCheck size={26} />}
//             step="STEP 4"
//             title="Clinical Insight"
//             desc="Results presented for clinical decision support and academic evaluation."
//           />

//         </div>
//       </section>

//       {/* ================= FEATURE HIGHLIGHTS ================= */}
//       <section className="bg-white py-16">
//         <div className="grid md:grid-cols-3 gap-16 text-center max-w-6xl mx-auto">

//           <div>
//             <h3 className="text-xl font-semibold text-teal-700 mb-2">
//               3 Disease Classes
//             </h3>
//             <p className="text-gray-500">
//               Sickle Cell Disease, Cystic Fibrosis,
//               Huntington's Disease
//             </p>
//           </div>

//           <div>
//             <h3 className="text-xl font-semibold text-teal-700 mb-2">
//               XGBoost Model
//             </h3>
//             <p className="text-gray-500">
//               High-accuracy gradient boosting for robust classification
//             </p>
//           </div>

//           <div>
//             <h3 className="text-xl font-semibold text-teal-700 mb-2">
//               SHAP Explainability
//             </h3>
//             <p className="text-gray-500">
//               Transparent, interpretable predictions for clinical trust
//             </p>
//           </div>

//         </div>
//       </section>

//       {/* ================= FOOTER ================= */}
//       <footer className="bg-gray-100 border-t border-gray-200 py-6 text-center text-gray-500 text-sm">
//         © 2026 Genetic Mutation Detection & Disease Prediction System —
//         Academic Research Project
//       </footer>
//     </div>
//   );
// }


// /* ================= STEP CARD ================= */
// function StepCard({ icon, step, title, desc }) {
//   return (
//     <div className="bg-white border border-gray-200 rounded-2xl p-8 shadow-sm">
      
//       <div className="bg-teal-100 text-teal-700 w-12 h-12 flex items-center justify-center rounded-xl mb-5">
//         {icon}
//       </div>

//       <p className="text-xs text-gray-400 font-semibold tracking-widest mb-2">
//         {step}
//       </p>

//       <h3 className="text-lg font-semibold text-gray-800 mb-2">
//         {title}
//       </h3>

//       <p className="text-gray-500 text-sm leading-relaxed">
//         {desc}
//       </p>

//     </div>
//   );
// }


import React from "react";
import { Link } from "react-router-dom";
import {
  FlaskConical,
  Brain,
  BarChart3,
  ShieldCheck,
  ArrowRight,
} from "lucide-react";

export default function HomePage() {
  return (
    <div className="min-h-screen bg-gray-50">

      {/* ================= HERO SECTION ================= */}
      <section className="bg-linear-to-r from-teal-700 to-emerald-600 text-white text-center py-32 px-6">
        <div className="max-w-4xl mx-auto">

          <div className="inline-flex items-center gap-2 bg-white/20 px-6 py-2 rounded-full text-sm backdrop-blur-md mb-8">
            🔬 Machine Learning–Powered Genomic Analysis
          </div>

          <h1 className="text-5xl font-bold leading-tight mb-8">
            Genetic Mutation Detection &
            <br />
            Disease Prediction System
          </h1>

          <p className="text-lg text-gray-100 leading-relaxed mb-10 max-w-3xl mx-auto">
            Leveraging XGBoost and Explainable AI (SHAP) to detect genetic
            mutations and classify diseases — Sickle Cell Disease,
            Cystic Fibrosis, and Huntington's Disease — with
            transparent, interpretable predictions.
          </p>

          <div className="flex justify-center gap-6">
            
            <Link
              to="/predict"
              className="flex items-center gap-2 bg-emerald-400 text-white px-8 py-3 rounded-xl font-semibold hover:opacity-90 transition"
            >
              <FlaskConical size={18} />
              Start Prediction
            </Link>

           <Link
  to="/about"
  className="flex items-center gap-2 border border-white px-8 py-3 rounded-xl font-semibold hover:bg-white hover:text-teal-700 transition"
>
  Learn More
  <ArrowRight size={16} />
</Link>

          </div>
        </div>
      </section>

      {/* ================= HOW IT WORKS ================= */}
      <section className="bg-gray-100 py-24 px-12">
        <div className="text-center mb-16">
          <h2 className="text-3xl font-bold text-gray-800 mb-3">
            How It Works
          </h2>
          <p className="text-gray-500">
            A step-by-step overview of the prediction pipeline
          </p>
        </div>

        <div className="grid md:grid-cols-4 gap-10 max-w-7xl mx-auto">

          <StepCard
            icon={<FlaskConical size={26} />}
            step="STEP 1"
            title="Input Genetic Data"
            desc="Enter mutation features: risk score, mutation status, and clinical significance."
          />

          <StepCard
            icon={<Brain size={26} />}
            step="STEP 2"
            title="ML Prediction"
            desc="XGBoost model classifies the mutation into SCD, CF, or Huntington's Disease."
          />

          <StepCard
            icon={<BarChart3 size={26} />}
            step="STEP 3"
            title="SHAP Explanation"
            desc="Explainable AI reveals which features drove the prediction."
          />

          <StepCard
            icon={<ShieldCheck size={26} />}
            step="STEP 4"
            title="Clinical Insight"
            desc="Results presented for clinical decision support and academic evaluation."
          />

        </div>
      </section>

      {/* ================= FEATURE HIGHLIGHTS ================= */}
      <section className="bg-white py-16">
        <div className="grid md:grid-cols-3 gap-16 text-center max-w-6xl mx-auto">

          <div>
            <h3 className="text-xl font-semibold text-teal-700 mb-2">
              3 Disease Classes
            </h3>
            <p className="text-gray-500">
              Sickle Cell Disease, Cystic Fibrosis,
              Huntington's Disease
            </p>
          </div>

          <div>
            <h3 className="text-xl font-semibold text-teal-700 mb-2">
              XGBoost Model
            </h3>
            <p className="text-gray-500">
              High-accuracy gradient boosting for robust classification
            </p>
          </div>

          <div>
            <h3 className="text-xl font-semibold text-teal-700 mb-2">
              SHAP Explainability
            </h3>
            <p className="text-gray-500">
              Transparent, interpretable predictions for clinical trust
            </p>
          </div>

        </div>
      </section>

      {/* ================= FOOTER ================= */}
      <footer className="bg-gray-100 border-t border-gray-200 py-6 text-center text-gray-500 text-sm">
        © 2026 Genetic Mutation Detection & Disease Prediction System —
        Academic Research Project
      </footer>
    </div>
  );
}


/* ================= STEP CARD ================= */
function StepCard({ icon, step, title, desc }) {
  return (
    <div className="bg-white border border-gray-200 rounded-2xl p-8 shadow-sm">
      
      <div className="bg-teal-100 text-teal-700 w-12 h-12 flex items-center justify-center rounded-xl mb-5">
        {icon}
      </div>

      <p className="text-xs text-gray-400 font-semibold tracking-widest mb-2">
        {step}
      </p>

      <h3 className="text-lg font-semibold text-gray-800 mb-2">
        {title}
      </h3>

      <p className="text-gray-500 text-sm leading-relaxed">
        {desc}
      </p>

    </div>
  );
}