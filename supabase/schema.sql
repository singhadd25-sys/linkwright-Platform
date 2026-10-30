-- Reference schema for Linkwright researcher profiles.
-- The matching RLS schema has already been applied to Supabase project ppxrlrevcitafcksajpi.
create table if not exists public.researcher_profiles (
  user_id uuid primary key references auth.users(id) on delete cascade,
  display_name text,
  created_at timestamptz not null default now(),
  updated_at timestamptz not null default now()
);
alter table public.researcher_profiles enable row level security;
revoke all on table public.researcher_profiles from anon;
grant select, insert, update on table public.researcher_profiles to authenticated;
drop policy if exists "Researchers can view their own profile" on public.researcher_profiles;
create policy "Researchers can view their own profile" on public.researcher_profiles
  for select to authenticated using ((select auth.uid()) = user_id);
drop policy if exists "Researchers can create their own profile" on public.researcher_profiles;
create policy "Researchers can create their own profile" on public.researcher_profiles
  for insert to authenticated with check ((select auth.uid()) = user_id);
drop policy if exists "Researchers can update their own profile" on public.researcher_profiles;
create policy "Researchers can update their own profile" on public.researcher_profiles
  for update to authenticated using ((select auth.uid()) = user_id)
  with check ((select auth.uid()) = user_id);
