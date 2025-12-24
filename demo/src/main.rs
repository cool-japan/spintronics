//! Spintronics Interactive Web Demonstrations
//!
//! HTMX + Axum powered real-time physics simulations
//!
//! Copyright (c) 2025 COOLJAPAN OÜ (Team KitaSan)
//! Licensed under MIT OR Apache-2.0

use std::sync::Arc;

use axum::response::IntoResponse;
use axum::routing::{get, post};
use axum::Router;
use tower_http::cors::CorsLayer;
use tower_http::services::ServeDir;
use tower_http::trace::{DefaultMakeSpan, TraceLayer};
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;

mod demos;
use demos::*;

mod templates;

#[cfg(test)]
mod tests;

/// Application state
#[derive(Clone)]
struct AppState {
    // Add shared state here if needed
}

#[tokio::main]
async fn main() {
    // Initialize tracing
    tracing_subscriber::registry()
        .with(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| "spintronics_demo=debug,tower_http=debug,axum=trace".into()),
        )
        .with(tracing_subscriber::fmt::layer())
        .init();

    let state = Arc::new(AppState {});

    // Build router
    let app = Router::new()
        // Main pages
        .route("/", get(index))
        .route("/llg", get(llg_page))
        .route("/spin-pumping", get(spin_pumping_page))
        .route("/materials", get(materials_page))
        .route("/skyrmion", get(skyrmion_page))
        // HTMX API endpoints
        .route("/api/llg/simulate", post(llg_simulate))
        .route("/api/llg/step", post(llg_step))
        .route("/api/spin-pumping/calculate", post(spin_pumping_calculate))
        .route("/api/materials/compare", get(materials_compare))
        .route("/api/skyrmion/visualize", post(skyrmion_visualize))
        // Static files
        .nest_service("/static", ServeDir::new("demo/static"))
        .nest_service("/public", ServeDir::new("demo/public"))
        .with_state(state)
        .layer(CorsLayer::permissive())
        .layer(
            TraceLayer::new_for_http()
                .make_span_with(DefaultMakeSpan::default().include_headers(true)),
        );

    // Run server
    let listener = tokio::net::TcpListener::bind("127.0.0.1:3000")
        .await
        .unwrap();
    tracing::info!("🚀 Spintronics Demo Server listening on http://127.0.0.1:3000");
    axum::serve(listener, app).await.unwrap();
}

/// Home page
async fn index() -> impl IntoResponse {
    templates::IndexTemplate {}.into_response()
}

/// LLG dynamics page
async fn llg_page() -> impl IntoResponse {
    templates::LlgTemplate {}.into_response()
}

/// Spin pumping page
async fn spin_pumping_page() -> impl IntoResponse {
    templates::SpinPumpingTemplate {}.into_response()
}

/// Materials explorer page
async fn materials_page() -> impl IntoResponse {
    templates::MaterialsTemplate {}.into_response()
}

/// Skyrmion simulator page
async fn skyrmion_page() -> impl IntoResponse {
    templates::SkyrmionTemplate {}.into_response()
}
